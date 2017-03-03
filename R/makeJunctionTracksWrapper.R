makeJunctionTracksWrapper <- function(bamFiles, genome,
                                      cores=1L,
                                      seqlev=NULL,
                                      min_junction_count=2,
                                      bedDir=".",
                                      juncDir=".",
                                      chrom_sizefile,
                                      trackName) {
    ## this wrapper makes bed12 files for junction and convert bed files
    ## to bigBed files for UCSC genome brower
    ## the main directory of juncDir should be
    ## "/fh/fast/tapscott_s/pub/tapscott/ucsc/junctions" 
   

    #' (1) santity check
    bdFiles <- makeJunctionBed(bamFiles=bamFiles, cores=cores, seqlev=seqlev,
                               min_junction_count=min_junction_count,
                               genome=genome, verbose=TRUE,                   
                               ignore.strand=TRUE, outdir=bedDir)
    convertBedToBigBed(bedFiles, sample_name=NULL, cores=cores,
                       chrom_sizefile, outdir=juncDir)
    createHubTrackLine.junc.bed(bbDir=juncDir,
                                trackName=trackName,
                                pattern="\\.bb$")
    
}
                                      

