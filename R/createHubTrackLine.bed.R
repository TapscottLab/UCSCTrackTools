createHubTrackLine.junc.bed <- function(bbDir,
                                        trackName,
                                        HubFileName=NULL,
                                        pattern="\\.bb$") {

    type <- "bigBed 12"

    
    ## take space out of hubName
    #hubName <- sub(" ", "_", hubName, fixed=TRUE)
    trackName <- gsub("[[:space:]]", "-", trackName)
    projectName <- basename(bbDir)
    #' sanity track
    if (!file.exists(bbDir)) stop(bbDir, " does not exist.")
    
    if (is.null(HubFileName))
        HubFileName <- paste0(projectName, "_HubTrackLines.txt")

    url <- "http://tapscott:FSHD@xfiles.fhcrc.org:7007/ucsc/tapscott/junctions"
    flist.bed <- list.files(bbDir,
                           full.name=TRUE, pattern=pattern)

    url.bed <- file.path(url, projectName, basename(flist.bed))
 
    if (length(flist.bed) == 0) stop("BED files do not exist.")


    HubTracks <- function(name, url.bed, type) {
        cat("track", name,  "\n")
        cat("compositeTrack on\n")
        cat("shortLabel", name, "\n")
        cat("longLabel", name, "\n")
        cat("allButteronPair on\n")
        cat("itemRgb on\n")
        cat("useScore 1\n")
        cat("type", type, "\n")
        cat("dragAndDrop on\n")


        cat("\n")    
          
        for(i in 1:length(url.bed)) {
            trackname <- basename(url.bed[i])
            cat("track", trackname, "\n")
            cat("type", type ,"\n")
            cat("shortLabel", trackname, "\n");
            cat("longLabel", trackname, "junction", "\n")
            cat("parent", name, "\n")
            cat("bigDataUrl", url.bed[i], "\n")
            cat("\n")
        }   
    }

    tFile <- file.path(bbDir, HubFileName)
    message("Export ", tFile)
    file.create(tFile)
    sink(tFile)
    HubTracks(name=trackName, url.bed=url.bed, type=type)
    sink()
 
}
