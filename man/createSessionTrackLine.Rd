\name{createSessionTrackLine.BW}
\alias{createSessionTrackLine.BW}

\title{
Create UCSC Genome Browser Session Track Line
}

\description{
Create a session track line txt file 
}
\usage{
createSessionTrackLine.BW(projectName, sessionFileName = NULL,
                       type = "bigWig", sample.group,
		       pattern="\\.bw$", group.col = NULL,
		       brewer = "Set2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{projectName}{
  The sub-directory in /home/tapscott/ucsc/bigWig where the bigWig files are stored.
}
  \item{sessionFileName}{
  A character string indicating the name of the session track file. If NULL, a txt file
  named "sessionTrackLines.txt" is created and saved in the current directory.
}
  \item{type}{
  "bigWig".}
  
  \item{sample.group}{
  A factor indicating the sample grouping and the color of the track will be assigned
  accordingly.}
  
  \item{pattern}{An optional regular expression.  Only file names which match
          the regular expression will be returned. The default if the file names end
	  up with .bw.}
	  
  \item{group.col}{
  A character string indicating the color for the track. If "NULL", the color patette is 
  signed by the brewer arguement.}
  
  \item{brewer}{
  A character string indicating the name of palette. See RColorBrewer's brewer.pal
  function.}
  
}
\details{
  The bigWig files must be placed in /home/tapscott/ucsc/bigWig.
}

\value{

}

\author{
Chao-Jen Wong <cwon2@fhcrc.org>
}

\examples{

\dontrun{
setwd("/home/tapscott/ucsc/bigWig/hg19_Dux4_RD_inducible")
getwd()
createSessionTrackLine.BW(projectName="hg19_Dux4_RD_inducible", sample.group=c(1,2))
}

}

