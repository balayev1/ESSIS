############################################### Execute_fishHook.R
################ This R script contains a function to execute fishHook, a package that enables statistical analysis of coding and non-coding mutation recurrence through GLM
################ This script also allows for inclusion of covariates into fisHook analysis such as DNA replication time, dinucleotide and trinucleotide percentages and read depth across genomic 
################ regions. See Generate_covariates4fishHook.sh file if you need to create the covariate tables.

### install required packages
require("optparse")
require("fishHook")    
require("gTrack")
require("rtracklayer")
require("skitools")

### set the arguments
option_list = list(
    make_option(c("-x", "--tiles-file"), type="character", default=NULL,
                help="Full path to precomputed genome tiles GRanges file in RDS format. If not provided, tiles will be computed automatically from mutation seqinfo.", metavar="character"),
    make_option(c("-m", "--mutation-file"), type="character", default=NULL, 
                help="full path to mutation file in RDS format", metavar="character"),
    make_option(c("-e", "--eligible-regions"), type="character", default=NULL, 
                help="full path to file with eligible  regions in RDS format (see FishHook tutorial for the right format)", metavar="character"),
    make_option(c("-d", "--id-column"), type="character", default=NULL, 
                help="Sample ID column in your mutation file", metavar="character"),
    make_option(c("-c", "--num-cores"), type="numeric", default=1, 
                help="Number of cores to run fishHook", metavar="numeric"),
    make_option(c("-p", "--id-cap"), type="numeric", default=1e9, 
                help="Maximum number of mutations allowed per sample per region", metavar="numeric"),
    make_option(c("-r", "--replication-time"), type="character", default=NULL, 
                help="full path to file including replication times of genomic regions in RDS format. MUST BE GRanges object. See tutorial for format details.", metavar="character"),
    make_option(c("-n", "--dinuc"), type="character", default=NULL, 
                help="full path to file including dinucleotide percentages of genomic regions in RDS format. MUST BE GRanges object. See tutorial for format details.", metavar="character"),
    make_option(c("-t", "--trinuc"), type="character", default=NULL, 
                help="full path to file including trinucleotide percentages of genomic regions in RDS format. MUST BE GRanges object. See tutorial for format details.", metavar="character"),
    make_option(c("-v", "--read-depth"), type="character", default=NULL, 
                help="full path to file including average read depths of genomic regions in RDS format. MUST BE GRanges object. See tutorial for format details.", metavar="character"),
    make_option(c("-b", "--h3k27ac-depth"), type="character", default=NULL, 
                help="full path to file including average H3K27ac read depths of genomic regions in RDS format. MUST BE GRanges object. See tutorial for format details.", metavar="character"),
    make_option(c("-f", "--h3k4me3-depth"), type="character", default=NULL, 
                help="full path to file including average H3K4me3 read depths of genomic regions in RDS format. MUST BE GRanges object. See tutorial for format details.", metavar="character"),
    make_option(c("-i", "--h3k36me3-depth"), type="character", default=NULL, 
                help="full path to file including average H3K36me3 read depths of genomic regions in RDS format. MUST BE GRanges object. See tutorial for format details.", metavar="character"),
    make_option(c("-j", "--h3k27me3-depth"), type="character", default=NULL, 
                help="full path to file including average H3K27me3 read depths of genomic regions in RDS format. MUST BE GRanges object. See tutorial for format details.", metavar="character"),
    make_option(c("-s", "--tile-size"), type="numeric", default=5e3, 
                help="Split size of genome tiles ", metavar="numeric"),
    make_option(c("-o", "--output-file"), type="character", default=NULL, 
                help="full path to output file to save fishHook output", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

run_fishook <- function(mutation_file, eligible_regions_file, tiles_file, tile_size, id_col, id_cap, num_cores, repli_time, dinuc, trinuc, read_depth, h3k27ac_rd, h3k4me3_rd, h3k36me3_rd, h3k27me3_rd, output_file){
    if (is.null(mutation_file)){
        stop("Please provide mutation file in RDS format\n")
    } else {
        mutations <- readRDS(mutation_file)
        message("Loaded ", length(mutations), " mutations")
    }

    if (is.null(output_file)) stop("Output file path is required.")
    
    ### Load eligible regions file if provided
    if (!(is.null(eligible_regions_file))){
        eligible_regions_file <- readRDS(eligible_regions_file)
    }

    ### Set tiles file
    if (!is.null(tiles_file)) {
        if (!file.exists(tiles_file)) {
            stop("Tiles file provided does not exist.")
        }
        tiles <- readRDS(tiles_file)
    } else {
        tiles <- gr.tile(seqinfo(mutations), tile_size)
    }

    ### Generate FishHook object
    fish.run = Fish(hypotheses = tiles, 
    events = mutations, 
    eligible = eligible_regions_file,
    idcol = id_col, 
    use_local_mut_density = TRUE, 
    mc.cores = num_cores, 
    idcap = id_cap)


    ### Add covariates if provided
    load_covariate <- function(path, field = NULL, name = NULL) {
    if (is.null(path)) return(NULL)
    gr <- readRDS(path)
    if (!is.null(name)) {
        Cov(data = gr, name = name)
    } else {
        Cov(data = gr, field = names(mcols(gr)))
    }
    }

    cov1 <- load_covariate(repli_time, name = "ReplicationTiming")
    cov2 <- load_covariate(dinuc)
    cov3 <- load_covariate(trinuc)
    cov4 <- load_covariate(read_depth, name = "AveReadDepth")
    cov5 <- load_covariate(h3k27ac_rd)
    cov6 <- load_covariate(h3k4me3_rd)
    cov7 <- load_covariate(h3k36me3_rd)
    cov8 <- load_covariate(h3k27me3_rd)
    
    covs <- c(cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8)
    if (!is.null(covs)){
        fish.run$covariates <- covs
    }

    ### Run the model
    fish.run$score()

    ### Add observed - expected difference
    res <- as.data.frame(fish.run$res)
    res$obs_minus_exp <- floor(res$count - res$count.pred)
    res$obs_minus_exp <- pmax(res$obs_minus_exp, 0)

    ### Save the final fishHook mutation matrix
    write.table(res, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)


    ### Print summary of the model fit
    if (!is.null(fish.run$model)) {
        message("Model summary:\n", capture.output(summary(fish.run$model), collapse = "\n"))
    }

    ### Print number of regions with FDR < 0.05
    message("Number of regions with FDR < 0.05: ", sum(fish.run$res$fdr < 0.05))
}

run_fishook(tiles_file = opt$'tiles-file',
    mutation_file = opt$'mutation-file',
    eligible_regions_file = opt$'eligible-regions',
    tile_size = opt$'tile-size',
    id_col = opt$'id-column',
    id_cap = opt$'id-cap',
    num_cores = opt$'num-cores',
    repli_time = opt$'replication-time',
    dinuc = opt$'dinuc',
    trinuc = opt$'trinuc',
    read_depth = opt$'read-depth',
    h3k27ac_rd = opt$'h3k27ac-depth',
    h3k4me3_rd = opt$'h3k4me3-depth',
    h3k36me3_rd = opt$'h3k36me3-depth',
    h3k27me3_rd = opt$'h3k27me3-depth',
    output_file = opt$'output-file')
