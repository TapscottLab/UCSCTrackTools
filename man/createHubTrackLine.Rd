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
                              hubName=NULL,
                              HubFileName = NULL,
                              col = c(102, 194, 165),
                              pattern="\\.bw$")
			      
    createHubTrackLine.BW(projectName, hubName=NULL,
                          HubFileName=NULL,
                          col = c(102, 194, 165), 
                          pattern="\\.bw$")
			  
    createHubTrackLine.junc.bed(projectName, hubName=NULL,
                                HubFileName=NULL,
                                pattern="\\.bb$")
				
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{projectName}{
  The sub-directory in /home/tapscott/ucsc/bigWig where the bigWig files are stored.
}
  \item{hubName}{Parent track name appears on the Hub.}

  \item{HubFileName}{
  A character string indicating the name of the hub track file. If NULL, a txt file
  named "$projectName__HubTrackLines.txt" is created and saved in the current directory.
}
  
  \item{col}{Color default to RGB:c(102, 194, 165) }

  \item{pattern}{The pattern of the bigWig file extension.}
  
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
                                      hubName="Dux4_RD_inducible")
}

}

