\name{makeCoverageTracksWrapper}
\alias{makeCoverageTracksWrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Coverage Tracks Wrapper
}
\description{
A wrapper callingconvertBamToBigWig() to make bigWig (.bw) coverage files of
RNA-seq data and save them to a designated destination.
}
\usage{
makeCoverageTracksWrapper(bamFiles, bwDir, singleEnded = TRUE,
                    NH.weight = FALSE, cores = 1, trackName,
		    TrackFileName = NULL, col = c(102, 194, 165),
		    pattern = "\\.bw$", shortLabel = NULL,
		    longLabel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamFiles}{
    charactor strings of the bam files subject to be converted to coverage bigWig files.
}
  \item{bwDir}{
    a character string specifying the destination of the bigWig files to be saved.
    For Tapscott lab, it HAS to be "/fh/fast/tapscott_s/pub/tapscott/ucsc/bigWig/..."
}
  \item{singleEnded}{
    logical indicating whether the reads are single ended.
}
  \item{NH.weight}{
    logical specifying whether the coverage should be weighted by the number of reported
    multiple alignments (NH column in bam files).
}
  \item{cores}{
    a numerical indicating how many cores is available.
}
  \item{trackName}{
    a character string specifying the name of the track to be appeared on UCSC genome browser.
}
  \item{TrackFileName}{
    a character string specifying the track line txt file. Default to trackName_HubTrackLines.txt
}
  \item{col}{
    a vector of three numerics indicating the rgb color.
}
  \item{pattern}{
    a character string specifing the default to "\\.bw$". The patter
}
  \item{shortLabel}{
    a character string indicating the short label attribute of the track.
}
  \item{longLabel}{
    a character string indicating the long label attribute of the track 
}
}

\value{
No value returns. The .bw files and a track line txt file will be saved to a designated destination indicated by bwDir parameter.
}

\author{
Chao-Jen Wong <cwon2@fredhutch.org>
}



\seealso{
\code{\link{convertBamToBigWig}}
}
\examples{
## The following case is Ashlee's canine RNA-seq samples. 
\dontrun{
pkgDir <- "/fh/fast/tapscott_s/R/RNA-Seq/canFam3.DuxFamily"
dataDir <- file.path(pkgDir, "data")
ngsDir <- "/shared/ngs/illumina/atyler/170208_SN367_0857_AHF37NBCXY"
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern="\\.bam$", full.names=TRUE)

bamFiles <- bamFiles[c(1, 5, 6)]

ucscDir <- "/fh/fast/tapscott_s/pub/tapscott/ucsc/bigWig"
bwDir <- file.path(ucscDir, "canFam3_CnMb_luciferase")
trackName <- "CnMb_luciferase_RNAseq"

library(RColorBrewer)
col <- col2rgb(brewer.pal(4, "Set2"))[, 1]
makeCoverageTracksWrapper(bamFiles=bamFiles,
                          bwDir=bwDir, #only one character
                          cores=1L, trackName=trackName,
                          col=col)
			  
setwd(bwDir)
dir()
}

}

