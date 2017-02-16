
getItemRgb <- function(score, pal="blue") {
    ## use only eight color - green series
    i <- cut(score, c(0, 5, 10, 20, 30, 50, 100, 200, Inf))
    col <- RColorBrewer::brewer.pal(9, "Blues")[2:9]
    itemRgb <- col[i]      
}

makeJunctionBed <- function(bamFiles, cores=1, seqlev=NULL,
                            min_junction_count=2,
                            genome=NULL, verbose=TRUE,                   
                            ignore.strand=TRUE, outdir=".") {
    require(rtracklayer)
    require(Rsamtools)
    if (is.null(genome)) stop("Must assign genome, i.e., mm10")
    if (!file.exists(outdir)) stop(outdir, "does not exists.")
    
    ## bamFiles must be just charactor vector
    si <- seqinfo(Rsamtools::BamFileList(bamFiles))

    sl <- seqlevels(si)
    st <- rep(c("+", "-"), rep(length(si), 2))

    which <- GRanges(sl, IRanges(1, seqlengths(si)[sl]), strand="*")
    if (!is.null(seqlev)) {
        ## make sure seqlev is a subset of seqlevels(which)
        which <- keepSeqlevels(which, value=seqlev)
    }
    
    list_which <- split(which, seq_along(which))

    #' now only support single-Ended reads
    ignore.strand <- TRUE
    min_anchor <- 1
    ## min_junction_count <- 1
    paired_end <- FALSE

    for (bam_file in bamFiles) {
        if (verbose) message(bam_file)
        list_junctionsTrack <- mclapply(list_which, function(x) { #mclapply
            if (verbose)
                message(sprintf("%s:%s-%s(%s)",
                                seqnames(x), start(x), end(x), strand(x)))
            
            flag <- scanBamFlag(isSecondaryAlignment = FALSE)
            param <- ScanBamParam(flag = flag, tag = "XS", which = x)
            gap <- GenomicAlignments::readGAlignments(file = bam_file,
                                                  param = param)
        
            frag_exonic <- reduce(ranges(grglist(gap, drop.D.ranges = TRUE)))
            frag_intron <- ranges(GenomicAlignments::junctions(gap))
        
            ## predict junctions
            junctions <- unique(unlist(frag_intron)) + 1
            if (identical(0L, length(junctions))) {return()}
        
            score <- SGSeq:::junctionCompatible(junctions, frag_exonic,
                                        frag_intron, min_anchor)
            mcols(junctions) <- DataFrame(score=score)
            junctions <-
                junctions[which(mcols(junctions)$score >= min_junction_count)]
            if (identical(0L, length(junctions))) {return()}

            ## prepare for the track
            junctionsTrack <- GRangesForUCSCGenome(
                genome=genome,
                chrom=as.character(seqnames(x)),
                ranges=junctions,
                strand=strand(x))
        
            block_starts <- rep(1, length(junctionsTrack)*2)
            block_starts[seq.int(2, length(block_starts), by=2)] <-
                width(junctionsTrack)-1
            blocks <- IRanges(start=block_starts, width=2)
            blocks <- split(blocks, rep(1:length(junctionsTrack), each=2))
            itemRgb <- getItemRgb(mcols(junctions)$score)
        
            mcols(junctionsTrack) <- DataFrame(score=mcols(junctions)$score,
                                               thick=IRanges(
                                               start=start(junctionsTrack),
                                               end=end(junctionsTrack)),
                                               blocks=blocks,
                                               itemRgb=itemRgb) 
            junctionsTrack
            }, mc.cores=cores)

        names(list_junctionsTrack) <- NULL
        keep <- which(elementNROWS(list_junctionsTrack) > 0)
        list_junctionsTrack <- list_junctionsTrack[keep]
        sample_name <- sub(".bam", "", basename(bam_file))
        bed_file <- paste0(sample_name, ".bed")

        junctionsTrack <- do.call(c, list_junctionsTrack)
        junctionsTrack <- keepSeqlevels(junctionsTrack, value=seqlev)
        strand(junctionsTrack) <- "+"
        ## don't include trackLine because bigBed does not neet it
        if (verbose) message("Exporting ", file.path(outdir, bed_file))
        export(junctionsTrack,
               con=file.path(outdir, bed_file))
    }
    invisible()
}

## move to somewhere else
convertBedToBigBed <- function(bed_files, sample_name=NULL, cores=1,
                               chrom_sizefile, outdir=".") {
    require("BiocParallel")
    if (!file.exists(outdir))
        stop(outdir, "does not exist!")
    if (!all(file.exists(bed_files)))
        stop("Some bed files does not exist!")
        
    if (is.null(sample_name))
        sample_name <- sub(".bed", "", basename(bed_files))
    bb_files <- file.path(outdir, paste0(sample_name, ".bb"))
    curdir <- getwd()
    setwd("~/tapscott/bin")
    cmd <- sprintf("./bedToBigBed %s %s %s", bed_files, chrom_sizefile,
                   bb_files)
    bplapply(cmd, system, BPPARAM=MulticoreParam(worker=cores))
    setwd(curdir)
}
