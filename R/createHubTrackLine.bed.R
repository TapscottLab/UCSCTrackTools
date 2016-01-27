createHubTrackLine.junc.bed <- function(projectName, hubName=NULL,
                                  HubFileName=NULL,
                                  pattern="\\.bb$") {
    type <- "bigBed 12"
    if (is.null(hubName)) hubName=projectName
    
    if (!file.exists(file.path("/home/tapscott/ucsc/junctions", projectName)))
        stop(file.path("/home/tapscott/ucsc/junctions", projectName), " does not exist.")
    
    if (is.null(HubFileName))
        HubFileName <- file.path(getwd(), paste0(projectName, "_HubTrackLines.txt"))

    url <- "http://tapscott:FSHD@xfiles.fhcrc.org:7007/ucsc/tapscott/junctions"
    flist.bed <- list.files(file.path("/home/tapscott/ucsc/junctions",
                           projectName),
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
            trackname <- paste0(basename(url.bed[i], "_junc"))
            cat("track", trackname, "\n")
            cat("type", type ,"\n")
            cat("shortLabel", trackname, "\n");
            cat("longLabel", trackname, "junc", "\n")
            cat("parent", name, "\n")
            cat("bigDataUrl", url.bed[i], "\n")
           cat("\n")
        }   
    }

    tFile <- HubFileName
    file.create(tFile)
    sink(tFile)
    if (is.null(hubName)) hubName=projectName
    
    HubTracks(name=hubName, url.bed, type)
    sink()
 
}
