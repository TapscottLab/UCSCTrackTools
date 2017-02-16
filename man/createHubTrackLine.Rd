\name{createHubTrackLine}
\alias{createHubTrackLine.BW}
\alias{createHubTrackLine.junc.bed}
\alias{createSessionTrackLine.BW}

\title{
Create UCSC Genome Browser Hub Track Line Attributes for bigWig Files
}

\description{
Create a session track line txt file for bigWig files
}
\usage{
    createSessionTrackLine.BW(projectName,
                              trackName=NULL,
                              TrackFileName = NULL,
                              col = c(102, 194, 165),
                              pattern="\\.bw$",
			      shortLabel=NULL,
			      longLable=NULL)
			      
    createHubTrackLine.BW(projectName,
	                  trackName=NULL,
                          TrackFileName = NULL,
                          col = c(102, 194, 165),
                          pattern="\\.bw$",
			  shortLabel=NULL,
			  longLable=NULL)
			  
    createHubTrackLine.junc.bed(projectName,
                                trackName=NULL,
                                TrackFileName = NULL,
                                col = c(102, 194, 165),
                                pattern="\\.bw$",
			        shortLabel=NULL,
			        longLable=NULL)
				
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{projectName}{
  The sub-directory in /home/tapscott/ucsc/bigWig where the bigWig files are stored or /home/tapscott/ucsc/junctions where the bigBed files are stored.
}
  \item{trackName}{A character indicating the parent track attributes appears on the Hub.}

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
  The bigWig files must be placed in /home/tapscott/ucsc/bigWig.
}


\author{
Chao-Jen Wong <cwon2@fhcrc.org>
}

\examples{

\dontrun{
setwd("/home/tapscott/ucsc/bigWig/hg19_Dux4_RD_inducible")
getwd()
library(RColorBrewer)
col <- brewer.pal(3, "Set2")
col2rgb(col)
UCSCTrackTools::createHubTrackLine.BW(projectName="hg19_Dux4_RD_inducible",
                                      trackName="Dux4_RD_inducible",
				      longLabel="Dux4_RD_inducible RNA-Seq Coverage")
}

}

