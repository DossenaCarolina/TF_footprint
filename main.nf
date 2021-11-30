#!/usr/bin/env nextflow

/*
                              ==================================================
                               PIPELINE FOR FOOTPRINT ANALYSIS OF ATAC-SEQ DATA
                              ==================================================
Nextflow pipeline for the analysis of  transcription factor footprints in ATAC-seq datasets. Started in June 2020.
------------------------------------------------------------------------------------------------------------------
Pipeline overview:

1. Merge bam files from biological replicates
2. Shift reads in bam file
3. Create consensus bed file with all the peaks from the biological replicates
4. Obtain nucleotide sequences in fasta format within the peaks in the consensus
5. Search for sequences within the peaks that match the PWM of the TF selected
6. Determine if the putative TF binding sites are bound
------------------------------------------------------------------------------------------------------------------
*/


//Command to run 
//nextflow run main.nf --design [path/to/design/file.csv] {--tf [TF_NAME] | --tf_list [path/to/tf_list.txt]} -w [path/to/work/dir/]

//Validate inputs
if (!params.tf && !params.tf_list) {
println 'Any transcription factor name provided! Specify it using --tf (or --tf_list if a list file is provided) in the command line'
exit 1
}


//Read the file with the paths
Channel
    .fromPath(params.design)
    .splitCsv(header: true)
    .map { sample ->
    println([sample.name,sample.id,sample.bam])
    [sample.name,sample.id,file(sample.bam)] }
    .into { rep_to_group; 
	    in_bed_create }

in_bed_create
.collect { it[1] }
.set { in_consensus }

process CREATE_CONSENSUS_BED {

        label 'process_low'

        publishDir "${params.outdir}/peaksConsensus", mode: 'copy'
        container '/storage/data/MAP/hpcapps/opt/core/pipelines/ingm-base.0.5.2.simg'

        input:
        val id from in_consensus

        output:
        path 'consensus_peakset.bed' into in_obtain_fasta_sequences

        script:
        """
        ${HOME}/singularity/envs/footprinting/bin/Rscript ${HOME}/map/til-epigenetics/atac-seq/analysis/footprint/bin/consensus_create.R \\
        -s ${id.join(',')} \\
        -b ${params.boolean_consensus_rds} \\
        -c ${params.consensus_peakset} \\
        -o ./ \\
        """
}


rep_to_group
    .groupTuple(by: 0)
    .map { it ->  [ it[0], it[1].flatten(), it[2].flatten()] }
    .set { bam_mrep }
//    .println()

process MERGE_REP_BAM {

        tag "$name"
	label 'process_high'

	publishDir "${params.outdir}/${name}/mergedReplicate", mode: 'copy'
	container '/storage/data/MAP/hpcapps/opt/core/pipelines/ingm-chip-seq.0.1.7b0.5.2.simg'

	input:
	tuple val(name), val(ids), path(bams)  from bam_mrep

	output:
	tuple val(name), path("${name}.sorted.bam"), path("${name}.sorted.bam.bai") into bam_to_shift

	script:
	bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()
	def avail_mem = 3
        avail_mem = task.memory.toGiga()
  
        """
	source activate chip
        picard -Xmx${avail_mem}g MergeSamFiles \\
            ${'INPUT='+bam_files.join(' INPUT=')} \\
            OUTPUT=${name}.sorted.bam \\
            SORT_ORDER=coordinate \\
            VALIDATION_STRINGENCY=LENIENT \\
            TMP_DIR=tmp \\
	    USE_THREADING=true
        samtools index ${name}.sorted.bam
	"""
}


process SHIFT_BAM {
	
	tag "$name"
        label 'process_high'

	publishDir "${params.outdir}/${name}/shiftedBam", mode: 'copy'
	container '/storage/data/MAP/hpcapps/opt/core/pipelines/ingm-chip-seq.0.1.7b0.5.2.simg'

	input:
	tuple val(name), file(bam), file(bai) from bam_to_shift
	
	output:
	tuple val(name), file("${name}.shifted.sorted.bam"), file("${name}.shifted.sorted.bam.bai") into in_centipede
	
	script:
	"""
	source activate deeptools
	export TMPDIR=\$PWD
	alignmentSieve -b ${bam} -o ${name}.shifted.tmp.bam --ATACshift -p max
	
	samtools sort -@ 4 -O bam -o ${name}.shifted.sorted.bam ${name}.shifted.tmp.bam
	samtools index ${name}.shifted.sorted.bam
	"""
}
 
process OBTAIN_SEQUENCES {

	label 'process_medium'

	publishDir "${params.outdir}/atacPeaksFasta", mode: 'copy'
	container '/storage/data/MAP/hpcapps/opt/core/pipelines/ingm-chip-seq.0.1.7b0.5.2.simg'

	input:
	path bed from in_obtain_fasta_sequences

	output:
	file 'atac_peaks.fa'  into in_sequence_pwm_match

	script:
	"""
	source activate chip
	bedtools getfasta -fi ${params.genomes['hg38'].fasta}  -bed $bed -fo atac_peaks.fa
	""" 
}

//ch_bed_fasta= Channel.value("${params.outdir}/atacPeaksFasta/atac_peaks.fa")
ch_bed_fasta = in_sequence_pwm_match.first()

if (params.tf_list) {
Channel
    .fromPath(params.tf_list)
    .splitCsv()
//    .subscribe { println it }
    .set { ch_tf_list }

Channel
    .fromPath(params.jaspar_motifs)
    .splitCsv(header: true)
    .map { tf ->
//   println([tf.id,tf.name,tf.path,tf.length])
    [tf.id,tf.name,file(tf.path),tf.length] }   
    .combine(ch_tf_list)
    .filter { it[1] == it[4] }
    .map {it -> it[0..-2] }
    .set { in_fimo } 
} else {
Channel
    .fromPath(params.jaspar_motifs)
    .splitCsv(header: true)
    .map { tf ->
//   println([tf.id,tf.name,tf.path,tf.length])
    [tf.id,tf.name,file(tf.path),tf.length] }
    .filter { it[1] ==~ /${params.tf}/ }
    .set { in_fimo }
}

process MOTIF_MATCH {
        
        tag "$tf_name"
	label 'process_medium'

//	publishDir "${params.outdir}/${tf_name}_${id}/pwmMatchedSequences", mode: 'copy'
	container '/storage/data/MAP/hpcapps/opt/core/pipelines/ingm-base.0.5.2.simg'
	
	input:
	tuple val(id),val(tf_name),file(path),val(length) from in_fimo
	file peaks_fasta from ch_bed_fasta	

	output:
	tuple val(id),val(tf_name),path("${tf_name}_${id}.fimo.txt.gz"),val(length) into sites	

	
	script:
	"""
	${HOME}/singularity/envs/footprinting/bin/fimo --text --parse-genomic-coord ${path} ${peaks_fasta} | gzip > ${tf_name}_${id}.fimo.txt.gz
	"""
}

sites
.combine(in_centipede) 
.set { ch_in_centipede }

process RUN_CENTIPEDE {

        tag "${tf_name}:${name}"
	label 'process_medium'

	publishDir "${params.outdir}/${name}/${tf_name}_${id}", mode: 'copy'
	container '/storage/data/MAP/hpcapps/opt/core/pipelines/ingm-base.0.5.2.simg'

	input:
	tuple val(id),val(tf_name),path(fimo),val(length), val(name), file(bam), file(bai) from ch_in_centipede

	output:
	tuple file("*.txt"), file("*.pdf"), file("*.bed") into results

	script:
	"""
	${HOME}/singularity/envs/footprinting/bin/Rscript ${HOME}/map/til-epigenetics/atac-seq/analysis/footprint/bin/run_centipede.R \\
	-b ${bam} \\
	-f ${fimo} \\
	-o ./ \\
	-x ${name}_${tf_name} \\
	-l ${length} \\
	&>> ./centipede_Output_File.txt

	awk '\$11=="1"' ${name}_${tf_name}.postPrFitCentipede.txt | awk -v FS='\t' -v OFS='\t' '{print \$2, \$3, \$4, \$5, \$6}' > ${name}_${tf_name}.predictedregions.bed
	"""	
}
