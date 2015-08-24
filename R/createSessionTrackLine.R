#'
#' This script create session track lines
#' parameters: file names with address

createSessionTrackLine <- function(projectName, sessionFile=NULL, type="bigWig",
                                    sample.group, group.col=NULL, brewer="Set2") {
    if (is.null(sessionFile))
        sessionFile <- file.path(getwd(), "sessionTrackLines.txt")

    url <- "http://tapscott:FSHD@xfiles.fhcrc.org:7007/ucsc/tapscott/bigWig"
    flist.bw <- list.files(file.path("/home/tapscott/ucsc/bigWig/", projectName),
                           full.name=TRUE)
    url.bw <- file.path(url, projectName, basename(flist.bw))

    #require(dichromat)
    #require(biovizBase)
    if (!is.factor(sample.group)) sample.group <- factor(sample.group)

    mycol <- RColorBrewer::brewer.pal(length(levels(sample.group)), brewer)
    mycol <- col2rgb(mycol)
    group.col <- mycol[, sample.group]

    if (!is.null(group.col))
        group.col <- col2rgb(group.col)
    
    tFile <- sessionFileName
    file.create(tFile)
    sink(tFile)
    for (i in 1:length(flist.bw)) {
        trackname <- sub(".bw", "", basename(flist.bw[i]))
        trackLine <- new("GraphTrackLine",
                     type="bigWig",
                     name=trackname,
                     description=trackname,
                     visibility="full",
                     color=as.integer(group.col[, i]))

        toTrackLine <- paste0(as(trackLine, "character"), " ", "bigDataUrl=", url.bw[i])
        toTrackLine <- sub("bedGraph", "bigWig", toTrackLine)
        cat(toTrackLine)
        cat("\n\n")
    }
    sink()
}
