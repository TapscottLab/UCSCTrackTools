createHubTrackDbLines <-
    function(bwDir, trackName, trackDbFileName, bw.only=TRUE,
             nh.bw.only=FALSE, type="bigWig", color=c(204, 255, 255),
             dataUrl="http://tapscott:FSHD@xfiles.fhcrc.org:7007/ucsc/tapscott/bigWig",
             subLabels="RNAseq") {
    ## create two objects trackParent and subTracks
    ## Then convert to 
    if (!file.exists(bwDir)) stop(bwDir, "does not exist!")

    #trackname <- basename(bwDir)
    col <- paste(as.character(color), collapse=",")
    
    trackline <- sink(trackDbFileName)
    cat("track ", trackName)
    cat("\n")
    cat("compositeTrack on \n")
    cat("shortLabel ", trackName)
    cat("\n")
    cat("longLabel ", trackName)
    cat("\n")
    cat("allButteronPair on \n")
    cat("dragAndDrop on \n")
    cat("type bigWig 0 1.0 \n")
    cat("color ", col)
    cat("\n")
    cat("alwaysZero on \n\n")

    ## children track ".bw" only for now
    if (bw.only) {
        files <-list.files(bwDir)
        nh <- grep(".nh.bw", files)
        bwfiles <- files[-nh]
    }

    labels <- sub(".bw", "", bwfiles)
    tracks <- paste0(labels, "_", subLabels)
    
    for (i in 1:length(bwfiles)) {
        cat("track ", tracks[i])
        cat("\n")
        cat("type bigWig")
        cat("\n")
        cat("shortLabel ", labels[i])
        cat("\n")
        cat("longLabel ", labels[i])
        cat("\n")
        cat("parent", trackName)
        cat("\n")
        url <- file.path(dataUrl, trackName, bwfiles[i])
        cat("bigDataUrl ", url)
        cat("\n\n")
    }
    sink()

    message(trackDbFileName, " updated")
    invisible()
}

## examples
bwDir <- "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg19.MB135/inst/bw"
trackDbFileName <-
    "/fh/fast/tapscott_s/CompBio/RNA-Seq/hg19.MB135/inst/track/HubTrackDb.txt"
col <- col2rgb("#66C2A5")
createHubTrackDbLines(bwDir=bwDir, trackName="hg19_MB135_RNASeq",
                      trackDbFileName=trackDbFileName, color=col)


