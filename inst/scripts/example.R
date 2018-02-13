## RColorBrewer;  require(Rsamtools)
##    require(BiocParallel)
##    require(rtracklayer)

library(UCSCTrackTools)
ngsDir <- "/shared/ngs/illumina/acampbel/150408_SN367_0511_BHJ33TADXX"
bamDir <- file.path(ngsDir, "tophat")

tmp <- list.files(file.path(bamDir, "bam"), pattern="\\.bam$",
                  all.files=FALSE,
                  include.dirs=TRUE)
bamFiles <- file.path(bamDir, "bam", tmp)[1]

#'
#' Step 1: convert Bam files to bw files
#' 
res <- UCSCTrackTools::convertBamToBigWig(bamFiles, NH.weight=FALSE, output.dir=NULL)
res <- UCSCTRackTools::convertBamToBigWig(bamFiles, NH.weight=TRUE,
                                          output.dir=NULL)
convertBamToBigWig("dsakl.txt")

#'
#' Step 2: move the bamFiles to /home/tapscott/ucsc/bigWig/your file folder
#'

#' Step 3:
#' create session track lines: url, where is the bwfiles, where to
#' store the session track lines file
#' 
pkgDir <- "~/tapscott/IntegratedProject/UCSCTrackTools"
type <- "bigWig"
projectName <- "hg19_MB135_RNASeq"
sessionFileName <- file.path(pkgDir, "inst", "extdata", "sessionTRackLines.txt")

sample.group=c(rep("group1", 6), rep("group2", 6))

## problem: how to assign color to the group?
createSessionTrackLine(projectName, sessionFileName=NULL,
                        sample.group=sample.group, brewer="Set2")                           
