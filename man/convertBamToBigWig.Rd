\name{convertBamToBigWig}
\alias{convertBamToBigWig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert Bam files to BigWig Files using coverage on genome positoins.
}

\usage{
convertBamToBigWig(bamFiles, NH.weight = FALSE, output.dir = NULL)
}

\arguments{
  \item{bamFiles}{
  Aharactor strings the bam files.
}
  \item{NH.weight}{
  Logical indicating whether to weight the coverage by the NH tab.
}
  \item{output.dir}{
  Character string indicating the location of the output bigWig files.
}
}


\value{
List of coverage.
}

\author{
Chao-Jen Wong <cwon2@fhcrc.org>
}

\examples{
\dontrun{
ngsDir <- "/shared/ngs/illumina/acampbel/150408_SN367_0511_BHJ33TADXX"
bamDir <- file.path(ngsDir, "tophat")

tmp <- list.files(file.path(bamDir, "bam"), pattern="\\.bam$",
                  all.files=FALSE,
                  include.dirs=TRUE)
bamFiles <- file.path(bamDir, "bam", tmp)[1]

## convert Bam files to bw files
res <- UCSCTrackTools::convertBamToBigWig(bamFiles, NH.weight=FALSE)

}
}