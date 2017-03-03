\name{createHubTrackLine}
\alias{createHubTrackLine.BW}
\alias{createHubTrackLine.junc.bed}
\alias{createSessionTrackLine.BW}
\alias{ createHubTrackLine}
\alias{ createHubTrackLine.junc.bed}

\title{
Create UCSC Genome Browser Hub Track Line Attributes for bigWig Files
}

\description{
Create a hub or session track line txt file for bigWig files. The txt file is stored at current location.
}
\usage{
    createSessionTrackLine.BW(projectName,
                              trackName=NULL,
                              TrackFileName = NULL,
                              col = c(102, 194, 165),
                              pattern=".bw$",
			      shortLabel=NULL,
			      longLabel=NULL)
			      
    createHubTrackLine.BW(bwDir,
	                  trackName,
                          TrackFileName = NULL,
                          col = c(102, 194, 165),
                          pattern=".bw$",
			  shortLabel=NULL,
			  longLabel=NULL)
			  
    createHubTrackLine.junc.bed(projectName,
                                trackName=NULL,
                                TrackFileName = NULL,
                                col = c(102, 194, 165),
                                pattern=".bd$",
			        shortLabel=NULL,
			        longLabel=NULL)
				
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bwDir}{The location of the bigWig files. At Tapscott Lab, the bigWig file (.bw)
  should be placed at `/fh/fast/tapscott_s/pub/tapscott/ucsc/bigWig/your_project`.
  The junction bigBed files (.bb) should be placed at
  `/fh/fast/tapscott_s/pub/tapscott/ucsc/junctions/your_project`. }
  \item{projectName}{ 
  The sub-directory in /home/tapscott/ucsc/bigWig where the bigWig files are stored or
  /home/tapscott/ucsc/junctions where the bigBed files are stored.
}
  \item{trackName}{A character indicating the parent track attributes appears on the Hub.
  It should not have empty space.}

  \item{TrackFileName}{
  A character string indicating the name of the hub track file. If NULL, a txt file
  named "$projectName__HubTrackLines.txt" is created and saved in the current directory.
}
  
  \item{col}{Color default to RGB:c(102, 194, 165) }

  \item{pattern}{A character indicating the pattern of the bigWig file extension.}

  \item{shortLabel}{A character indicating the shortLabel attribute for the track lines.}

  \item{longLabel}{A character indicating the longLabel attribute for the track lines.}
  
}
\details{
  The bigWig files must be placed at `/fh/fast/tapscott_s/pub/tapscott/ucsc/bigWig/your_project`.
}


\author{
Chao-Jen Wong <cwon2@fhcrc.org>
}

\examples{
## this ucscDir is where the Tapscott Lab stores the bigWig files for the Hub
ucscDir <-  "/fh/fast/tapscott_s/pub/tapscott/ucsc/bigWig"
## define the sub-directory
bwDir <- file.path(ucscDir, "hg19_Dux4_RD_inducible")
library(RColorBrewer)
col <- col2rgb(brewer.pal(3, "Set2"))
createHubTrackLine.BW(bwDir=bwDir,
                      trackName="Dux4_RD_inducible",
     		      col=col[, 2],
		      longLabel="Dux4_RD_inducible RNA-Seq Coverage")

}

