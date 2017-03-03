createHubTrackLine.BW <- function(bwDir, trackName,
                                  TrackFileName=NULL,
                                  col = c(102, 194, 165), 
                                  pattern="\\.bw$",
                                  shortLabel=NULL,
                                  longLabel=NULL) {

    #' projectName: the sub-directory name used to store the bw files on the protal
    type <- "bigWig"
    projectName <- basename(bwDir)

    ## sanitize trackName, shortLabel and longLabel
    if (!file.exists(bwDir)) stop(bwDir, " does not exists.")
    trackName <- gsub("[[:space:]]", "-", trackName) ## does not allow empty space
    if (is.null(shortLabel)) shortLabel <- trackName
    if (is.null(longLabel)) longLabel <- trackName
    message("trackName: ", trackName)
    message("shortLabel: ", shortLabel)
    message("longLabel: ", longLabel)

    if (!file.exists(bwDir))
        stop(bwDir, " does not exist.")
    
    if (is.null(TrackFileName))
        TrackFileName <- file.path(getwd(),
                                   paste0(projectName, "_HubTrackLines.txt"))

    url <- "http://tapscott:FSHD@xfiles.fhcrc.org:7007/ucsc/tapscott/bigWig"
    flist.bw <- list.files(bwDir,
                           full.name=TRUE, pattern=pattern)
    if (length(flist.bw) == 0) stop("There is no bigWig files in ", bwDir)
    
    url.bw <- file.path(url, projectName, basename(flist.bw))

    if (length(flist.bw) == 0) stop("BigWig files do not exist.")


    HubTracks <- function(name, col, url.bw, shortLabel, longLabel) {
        cat("\n")
        cat("track ", name,  "\n")
        cat("compositeTrack on\n")
        cat("shortLabel ", shortLabel, "\n")
        cat("longLabel ", longLabel, "\n")
        cat("allButteronPair on\n")
        cat("dragAndDrop on\n")
        cat("type bigWig 0 1.0\n")
        cat("alwaysZero on\n")
        cat("color ", col, "\n")
        cat("\n")    
          
        for(i in 1:length(url.bw)) {
            trackname <- sub(".bw", "", basename(url.bw[i]))
            cat("track ", trackname, "\n")
            cat("type ", "bigWig\n")
            cat("shortLabel ", trackname, "\n");
            cat("longLabel ", trackname, " coverage", "\n")
            cat("parent ",name, "\n")
            cat("bigDataUrl ", url.bw[i], "\n")
           cat("\n")
        }   
    }

    tFile <- TrackFileName
    message("Creating ", tFile)
    file.create(tFile)
    sink(tFile)
    
    HubTracks(name=trackName, col, url.bw, shortLabel, longLabel)
    sink()
 
}
