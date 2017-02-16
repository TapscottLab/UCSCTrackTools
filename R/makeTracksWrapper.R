makeCoverageTracksWrapper <- function(bamFiles, bwDir, singleEnded=TRUE,
                                      NH.weight=FALSE, cores=1,
                                      trackName,
                                      TrackFileName = NULL, 
                                      col = c(102, 194, 165), pattern = "\\.bw$",
                                      shortLabel = NULL, longLabel = NULL) {
    #' This function
    #' (1) get coverage of a bam files and convert it bigWig (bw) files
    #' (2) create the track line files. The resulting files are placed at the "bwDir".
    #' This function wraps the convertBamToBigWig and createHubTrackLine.BW function
    #' and place the results in the desired directory.
    
    #' check if bamFiles exists
    if (!all(file.exists(bamFiles)))
        stop("Some of the bam files do not exists")

    #' check if the bwDir (bigWig) exists
    if (!file.exists(bwDir))
        stop(bwDir, " does  not exist")

    #' (1) get coverage from bamFiles and convert to bigWig files
    convertBamToBigWig(bamFiles, NH.weight=NH.weight, output.dir=bwDir,
                       cores=cores, singleEnded=singleEnded)
    curDir <- getwd()
    setwd(bwDir)
    projectName <- basename(bwDir)
    ucscDir <- dirname(bwDir)
    #' (2) create hub track lines
    createHubTrackLine.BW(bwDir=ucscDir, projectName=projectName,
                          trackName = trackName, TrackFileName = TrackFileName,
                          col=col, pattern=pattern, shortLabel=shortLabel, longLabel=longLabel)
    setwd(curDir)
}

    
