#!/usr/bin/env Rscript

##LOAD LIBRARIES

library(optparse)
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(CENTIPEDE))
suppressPackageStartupMessages(library(CENTIPEDE.tutorial))

## PARSE COMMAND-LINE PARAMETERS

option_list <- list(make_option(c("-b", "--bam_file"), type="character", default=NULL, help="ATAC-seq shifted bam file", metavar="path"),
                    make_option(c("-f", "--fimo_file"), type="character", default=NULL, help="File .txt.gz generated with fimo containing the nucleotide sequences within the peaks provided in the bed file matching the PWM of the selected motif", metavar="path"),
                    make_option(c("-p", "--pvalue"), type="double", default='1e-4', help="p value, level of significance"),
		    make_option(c("-s", "--flank_size"), type="integer", default='100', help="Defines the size of the region flanking upstream and downstream the PWM match sites to beconsidered to count the read start positions"),
		    make_option(c("-o", "--outdir"), type="character", default='./', help="Output directory", metavar="path"),
                    make_option(c("-x", "--outprefix"), type="character", default='sample_tf', help="Output prefix", metavar="string"),
make_option(c("-l", "--length"), type="integer", default=NULL, help="length of the motif"))


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

BamFile <- opt$bam_file
FimoFile <- opt$fimo_file
pval <- opt$pvalue
FlankSize <- opt$flank_size

## COUNT READ START POSITIONS

cen <- centipede_data(
  bam_file = BamFile,
  fimo_file = FimoFile,
  pvalue = pval,
  flank_size = FlankSize
)
 
##regions is a dataframe with one row for each PWM region.
##mat is a matrix with read counts for each PWM region.
## WRITE THE REGIONS AND COUNT DATA ON A FILE

regions <- as.data.frame(cen$regions)
mat <- as.matrix(cen$mat)

#RegionsFile <- file.path(opt$outdir,paste(opt$outprefix,".regions.txt",sep=""))
#write.table(regions,file=RegionsFile,quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)

MatFile <- file.path(opt$outdir,paste(opt$outprefix,".mat.txt",sep=""))
write.table(mat,file=MatFile,quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)

## PLOT READ START SITE COUNTS
#PlotFile <- file.path(opt$outdir,paste(opt$outprefix,".ReadStartSites.pdf",sep=""))
#if (file.exists(PlotFile) == FALSE) {
#    pdf(file=PlotFile,onefile=TRUE,width=7,height=7)
#
#plot <- plot(MatFile, xlab = "Position", ylab = "Read Start Sites", type = "h",
#     col = rep(c("blue", "red"), each = (2*FlankSize+13))
#abline(v = c(FlankSize, (FlankSize+13), (13+(3*FlankSize)), (26+(3*FlankSize))) + 0.5, lty = 2)
#abline(v = (13+(2*FlankSize)) + 0.5)
#print(plot)
#
#dev.off()
#}

## compute the posterior probability that a TF is bound at each peak

fit <- fitCentipede(
  Xlist = list(ATAC = cen$mat),
  Y = as.matrix(data.frame(
    Intercept = rep(1, nrow(cen$mat))
  ))
)

FitOutput <- cbind(regions, PostPr = fit$PostPr)
FitFile <- file.path(opt$outdir,paste(opt$outprefix,".postPrFitCentipede.txt",sep=""))
write.table(FitOutput,file=FitFile,quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)

## PLOT RESULTS

PlotFile <- file.path(opt$outdir,paste(opt$outprefix,".plotResults.pdf",sep=""))
if (file.exists(PlotFile) == FALSE) {
    pdf(file=PlotFile,onefile=TRUE,width=7,height=7)


heatmap <- imageCutSitesCombined(cen$mat[fit$PostPr == 1,])
footprint <- plotProfile(fit$LambdaParList[[1]], Mlen = opt$length)

print(heatmap)
print(footprint)

dev.off()
}
