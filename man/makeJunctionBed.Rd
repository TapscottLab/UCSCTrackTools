\name{makeJunctionBed}
\alias{makeJunctionBed}
\title{
Make Junction Bed Files from BAM files
}
\description{
Calculate juntions from BAMs files 
}
\usage{
makeJunctionBed(bamFiles, cores = 1, seqlev = NULL, genome = NULL,
                verbose = TRUE, ignore.strand = TRUE, outdir = ".")
}

\arguments{
  \item{bamFiles}{
  A character vector identifying BAM files that are used to find junctions.
  }
  \item{cores}{
  An integer indicating number of cores available for parallel processing.
  }
  \item{seqlev}{
  A character vector indicating the chromosome of interest.
  }
  \item{genome}{
  A character string identifying a genome, usually one assiagned by UCSC,
  like "mm10".
  }
  \item{verbose}{
  If "TRUE", generate message indicating progress.
  }
  \item{ignore.strand}{
  Logical indicating whether to ignore strand. At this moment, it is set to "TRUE" internally.
  }
  \item{outdir}{
  A string indicating the output directory.
  }
}

\details{
This function imports reads from BAM files and the method to calculate junction coutns is adapted from the SGSeq packages. 
}

\author{
Chao-Jen Wong <cwon2@fredhutch.org>
}

\seealso{createHubTrackLine.junc.bed}

\examples{
 \dontrun{
bam_file <- "/shared/ngs/illumina/jwhiddon/150918_SN367_0553_AH7YLFBCXX/tophat/bam/C2C12_mDuxCA_36h_15B_DOX1.bam"
seqlev <- paste0("chr", c(1:19, "M", "X", "Y"))
makeJunctionBed(bamFiles=bam_file, cores=4, seqlev=seqlev,
                genome="mm10",
                ignore.strand=TRUE, verbose=TRUE)
		
## if you want to view from UCSC Hub track (not session), you need to convert
## the bed file to bigBed format 
size_file <- "/fh/fast/tapscott_s/CompBio/mm10/mm10.chrom.sizes"
bedDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/mm10.C2C12.mDux/inst/bed"
UCSCTrackTools:::convertBedToBigBed(file.path(bedDir,
                                              "C2C12_mDuxCA_36h_15B_DOX1.bed"),
                                    chrom_sizefile=size_file,
				    cores=1)
## make hub trackline for bigBed or bed (15) files

 }
}

