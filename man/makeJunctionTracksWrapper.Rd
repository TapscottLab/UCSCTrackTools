\name{makeJunctionTracksWrapper}
\alias{makeJunctionTracksWrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Junction Tracks Wrapper 
}
\description{
A wrapper calls up three major functions to create junction bed file, junction
tracks bigBed files, and generate a track line files specifying the attribute of the
junciton track.
}
\usage{
makeJunctionTracksWrapper(bamFiles, genome, cores = 1L,
                          seqlev = NULL, min_junction_count = 2,
			  bedDir = ".", juncDir = ".",
			  chrom_sizefile,
			  trackName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamFiles}{
character strings  of the bam files.
}
  \item{genome}{
a character string specifying the genome of the subject.
}
  \item{cores}{
an integer indicating the number of cores.
}
  \item{seqlev}{
character strings specifying the sequence levels desired to be included in the
junction files.
}
  \item{min_junction_count}{
an integer specifying the minimal number of junction read.
}

  \item{bedDir}{
a character string indicating the destination of the junction bed (.bed) files.
}
  \item{juncDir}{
a character string indicating the destination of the junction bigBed (.bb) files.
The preferred location is "/fh/fast/tapscott_s/pub/tapscott/ucsc/junctions/...".
}
  \item{chrom_sizefile}{
a file name of the chromozon size file pertained to the genome. Those files can be
found in "/fh/fast/tapscott_s/CompBio/hg38, /mm10, and /canFam3". 
}
  \item{trackName}{
a character string specifying the name of the track to be
          appeared on UCSC genome browser.
}

}
\details{
The wrapper calls three functions: makeJunctionBed(), convertBedToBigBed(),
and createHubTrackLine.junc.bed(). See example below.
}

\author{
Chao-Jen Wong <cwon2@fhcrc.org}
}


\seealso{
\code{\link{makeJunctionBed}}, \code{\link{createHubTrackLine.junc.bed}}
}
\examples{

\dontrun{
pkgDir <- "/fh/fast/tapscott_s/R/RNA-Seq/canFam3.DuxFamily"
load(file.path(pkgDir, "data", "canine.ens.SE.rda"))
ucscDir <- "/fh/fast/tapscott_s/pub/tapscott/ucsc/junctions"
ngsDir <- "/shared/ngs/illumina/atyler/170208_SN367_0857_AHF37NBCXY"

## define parameters
bamDir <- file.path(ngsDir, "tophat", "bam")
bamFiles <- list.files(bamDir, pattern="\.bam$",
                       full.names=TRUE)[c(1,5, 6)]
seqlev <- seqlevels(canine.ens.SE)
chrom_sizefile <- "/fh/fast/tapscott_s/CompBio/canFam3/canFam3.chrom.sizes"
bedDir <- file.path(pkgDir, "inst", "junction_bed")
## bigBedDir is the location where files for the Hub reside.
bigBedDir <- file.path(ucscDir, "canFam3_CnMb_luciferase")
trackName <- "CnMb_luciferase_junc"

makeJunctionTracksWrapper(bamFiles, genome="canFam3",
                          cores=1L, seqlev=seqlev,
                          min_junction_count=2,
                          chrom_sizefile=chrom_sizefile,
                          bedDir=bedDir,
                          bigBedDir=bigBedDir,
                          trackName=trackName)
}
}

