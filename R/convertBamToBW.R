#'
#' This function converts BAM files to BigWig files by using the coverage on the
#' genome position. Return coverage.
#'

convertBamToBigWig <- function(bamFiles, NH.weight=FALSE, output.dir=NULL,
                               cores=1L, singleEnded=TRUE) {
    require(parallel)
    require(GenomicAlignments)
    require(rtracklayer)
    
    #NH.weight <- FALSE # not supporting NH weighted coverage
    if (is.null(output.dir)) 
        output.dir <- getwd()

    ## check if the directory exists
    if (!file.exists(output.dir)) stop("Output directory does not exists.")
    if (!all(file.exists(bamFiles))) stop("Bam files do not exists.")

    ## check the bam file is a bam file
    
    prefix <- sub(".bam", "", basename(bamFiles))
    param <- Rsamtools::ScanBamParam(tag="NH")

    res <- mclapply(bamFiles, function(filename, singleEnded, NH.weight=NH.weight) {
        print(message("Computing coverage for ", filename))
        
        if (singleEnded)
            reads <- readGAlignments(filename, param=param)
        else
            reads <- readGAlignmentPairs(filename, param=param)
        
            seqlevelsStyle(reads) <- "UCSC"
            cov <- GenomicAlignments::coverage(reads)

         if (NH.weight)  {
            weight <- 1.0/mcols(reads)$NH
            cov.nh <- coverage(reads, weight=weight)
            list(cov=cov, cov.nh=cov.nh, weight=weight, total=length(reads))
        }

        if (!NH.weight) {
            list(cov=cov, total=length(reads))
        }
        
    }, singleEnded=singleEnded, NH.weight=NH.weight,
                    mc.cores=cores, mc.preschedule=FALSE)

    names(res) <- prefix

    cov <- sapply(res, function(x) x$cov)
    names(cov) <- prefix

    ## save coverage as bw files
    mclapply(names(cov), function(x) {
        output <- file.path(output.dir, paste0(x,".bw"))
        print(message("Exporting BigWig file: ", output))
        rtracklayer::export(cov[[x]], con=output, format="BigWig")
    }, mc.cores=cores, mc.preschedule=FALSE)
    
    if (NH.weight) {
        cov.nh <- sapply(res, function(x) x$cov.nh)
        names(cov.nh) <- prefix
         ## save coverage with NH adjustment as bw files
        mclapply(names(cov.nh), function(x) {
            output <- file.path(putput, paste0(x, ".nh.bw"))
            print(message("Exporting BigWig file: ", ouptput))
            rtracklayer::export(cov.nh[[x]], output)
        }, mc.cores=cores, mc.preshcedule=FALSE)

    }

    return(res)

}
