# TF_footprint

# Footprint motif-centred analysis using CENTIPEDE

This pipeline is built using [Nextflow](https://www.nextflow.io) and can be used to perform computational footprint analysis on  ATAC-seq data. 

[[_TOC_]]

## Pipeline summary

1. Merge bam files from biological replicates ([`picard`](https://broadinstitute.github.io/picard/))
2. [Shift reads](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1642-2#Sec2) in bam file ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html))
3. Create consensus bed file with all the peaks from the samples in the design file ([`consensus_create.R`]())
4. Obtain nucleotide sequences in fasta format within the peaks in the consensus ([`BEDTools`](https://github.com/arq5x/bedtools2/))
5. Search for sequences within the peaks that match the PWM of the TF selected ([`FIMO`](http://meme-suite.org/doc/fimo.html))
6. Determine if the putative TF binding sites are bound ([`CENTIPEDE`](https://rdrr.io/rforge/CENTIPEDE/))

## Pipeline configuration

To run this pipeline you have to [install Nextflow](https://mapmachine.bioinfo.ifom.eu/gitlab/it/bioinformatics/-/wikis/home#installing-nextflow).

The pipelines works with [Singularity](https://sylabs.io/guides/3.5/user-guide/).

The execution of the following steps allows to recreate the conda environment in Singularity for high-reproducible results with CENTIPEDE: 


## Files required to run the pipeline


* A design file .csv with this structure: `biological_replicate_name,sample_id,/path/to/file.bam`
* A file .csv obtained from the JASPAR motif database with this structure: `motif_id,tf_name,/path/to/motif/file.meme,length_of_the_motif` for all the transcription factors presen in the Jaspar database 
* A consensus file in .bed format with all the peaks of all the samples of the dataset, produced as an output of the [nf-core atac-seq pipeline](https://mapmachine.bioinfo.ifom.eu/gitlab/map/pipelines/-/tree/master/atac-seq) (path to the file specified in `nextflow.config`)
* A sample boolean file in .rds format for the consensus peakset above mentioned (path to the file specified in `nextflow.config`). 


All the individual and non-redundant PFMs of motifs have been downloaded from the JASPAR CORE 2020 database:

     wget http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.zip

The csv file containing data and path to each motif file has been created using a  bash script that iterates over all of the files retrieving the necessary information. Here the script: 

    #!/bin/bash
    
    
    echo "id,name,path,length" >> jaspar_db.csv
    for i in JASPAR2020/*.meme
    do
    export name1=$(awk -vvar=2 -F" " 'NR==10 {print $var}' $i)
    export name2=$(awk -vvar=3 -F" " 'NR==10 {print $var}' $i)
    export length=$(awk -vvar=6 -F" " 'NR==11 {print $var}' $i)
    echo $name1,$name2,$i,$length >> jaspar_db.csv
    done
  

The pipeline comes with two bundled R scripts stored in `/bin`:
1. consensus_create.R
2. run_centipede.R

## Running the pipeline

The command for running the pipeline is as follows: 

    export TMPDIR='/path/to/tmp' 
    nextflow run main.nf --tf [TF_NAME] -w path/to/work/dir

**Required parameters**

* It is possible to specify either a single TF name/multiple TFs or a text file with a list of TFs in the command line. The analysis will be performed on the motfis of the TF specified. 
For the first option, specify `--tf [TF_NAME]` in the command line (e.g. `--tf 'CTCF'`). You can also specify more than one TF name at the same time and the pipeline will work as well; this is the correct syntax: `--tf "(CTCF|FOXP3)"`.
For the second option, specify `--tf_list [path/to/file.txt]` in the command line. It's foundamental to use either one or the other parameter altrernatively to run correctly the pipeline.

* **--design**: with this parameter you specify the path and name of your design file.
