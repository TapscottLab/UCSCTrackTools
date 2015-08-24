ngsDir <- "/shared/ngs/illumina/acampbel/150408_SN367_0511_BHJ33TADXX"
bamDir <- file.path(ngsDir, "tophat")

tmp <- list.files(file.path(bamDir, "bam"), pattern="\\.bam$",
                  all.files=FALSE,
                  include.dirs=TRUE)
bamFiles <- file.path(bamDir, "bam", tmp)[1]

## convert Bam files to bw files
res <- convertBamToBigWig(bamFiles, NH.weight=FALSE)
res <- convertBamToBigWig(bamFiles, NH.weight=TRUE)
convertBamToBigWig("dsakl.txt")
## move the bamFiles to /home/tapscott/ucsc/bigWig/your file folder

## create session track lines: url, where is the bwfiles, where to
## store the session track lines file

type <- "bigWig"
projectName <- "hg19_MB135_RNASeq"
sessionFileName <- "sessionTRackLines.txt"
sample.group=c(rep("group1", 6), rep("group2", 6))

createSessionTrackLines(projectName, sessionFile=NULL,
                        sample.group, brewer="Set2")                           
