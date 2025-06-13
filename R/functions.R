############################################################################
## Agshin Balayev
## University of Minnesota Neurosurgery  balay011@umn.edu
###############################################################################


#' @title parse_genomic_elements
#' @description
#'
#' Takes input of full path to GRanges of genomic elements of protein-coding genes (e.g. CDS, 5UTR, 3UTR, introns, promoters),
#' non-coding RNA genes (e.g. lncRNA exons/introns, miRNA, tRNA) and gene regulatory regions (e.g. superenhancers, enhancers, CTCF-binding motifs)
#' in RDS format. File is parsed based on the unique elements within the 'type' column.
#' Parsed files are later used as tiles in fishHook to identify excess mutations within each genomic element. Optional GTF provides intergenic regions
#' that excludes all regions in the input and VDJ sites of immunoglobulin and T-cell receptor genes. GTF must correspond to GRCh38 genome with chr1-22,X,Y
#' being considered.
#' Input RDS must have a metadata column named `"type"`.
#'
#' @param ge_file   Full path to GRanges RDS with genomic elements (required)
#' @param output_dir    Directory where parsed RDS files will be written (required) 
#' @param gtf_path  Full path to a GTF annotation to acquire intergenic regions (optional)
#' @return  Invisibly writes the parsed genomic elements within [output_dir]
#' @export
parse_genomic_elements <- function(ge_file, output_dir = NULL, gtf_path = NULL){

    # check for dependencies
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop("Please install the 'GenomicRanges' package")
    }
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        stop("Please install the 'rtracklayer' package")
    }

    # read input RDS file if file is provided and exists
    if (is.null(ge_file) || !file.exists(ge_file)){
        stop("Input RDS file missing")
    } else {
        ge.input <- readRDS(ge_file)
    }

    # split by unique values of 'type' column if exists
    if (!"type" %in% names(mcols(ge.input))) {
        stop("Input GRanges must have a 'type' metadata column")
    }
    types  <- unique(mcols(ge.input)$type)
    ge.list <- setNames(lapply(types, function(t) ge.input[mcols(ge.input)$type == t]), types)

    # import GTF file if provided and exists
    if (!is.null(gtf_path) && file.exists(gtf_path)){
        gtf <- import(gtf_path)

        # select VDJ sites of immunoglobulin and T-cell receptor genes
        exclude.types <- c("IG_V_gene","IG_J_gene","IG_D_gene","IG_V_pseudogene","IG_C_gene",
            "IG_J_pseudogene","IG_C_pseudogene","IG_pseudogene",
            "TR_C_gene","TR_J_gene","TR_V_gene","TR_V_pseudogene",
            "TR_D_gene","TR_J_pseudogene")
        excluded <- gtf[gtf$gene_type %in% exclude.types]

        # specify GRanges for GRCh38 chr1–22,X,Y
        chr.lengths <- c(248956422,242193529,198295559,190214555,181538259,170805979,
            159345973,145138636,138394717,133797422,135086622,133275309,
            114364328,107043718,101991189,90338345,83257441,80373285,
            58617616,64444167,46709983,50818468,156040895,57227415)
        seqs <- paste0("chr", c(1:22, "X", "Y"))
        genome <- GRanges(seqnames = seqs, ranges = IRanges(start=1, end=chr.lengths))
        seqinfo(genome) <- Seqinfo(seqs, chr.lengths, isCircular=rep(FALSE, length(chr.lengths)))

        # exclude all regions in the input RDS and VDJ sites of immunoglobulin and T-cell receptor genes
        subtract.regions <- reduce(c(excluded, ge.input))
        overlaps <- findOverlaps(genome, subtract.regions)
        grl <- extractList(subtract.regions, as(overlaps, "List"))
        intergenic <- do.call(c, psetdiff(genome, grl))

        # intergenic regions
        intergenic$type <- "INTERGENIC"
        ge.list[["INTERGENIC"]] <- intergenic
    }

    # save each element as RDS
    if (is.null(output_dir)) {
        stop("Output directory required to store parsed genomic elements")
    }

    if (!is.null(output_dir) && !dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    for (nm in names(ge.list)) {
        outfile <- file.path(output_dir, paste0(tolower(nm), ".Rds"))
        saveRDS(ge.list[[nm]], file = outfile)
    }

    invisible(NULL)

}


#' @title run_fishhook
#' @description
#'
#' Runs fishHook using each GRanges element set as hypotheses and each mutation file representative of each tumor cohort x mutation class as events.
#' Mutation classes are single-nucleotide variant (SNV) and insertion/deletion (INDEL). Mutation files must be included in text file with format:
#' column I: tumor cohort (i.e.tumor type), column II: mutation class, column III: full path to RDS GRanges mutations, column IV: path to RDS GRanges 
#' average read depth (optional). Covariates such as replication timing, dinucleotide percentages, trinucleotide percentages, average read depth 
#' and local mutation density can be included within fishHook model. 
#' All covariate RDS GRanges must have a single metadata column.
#' Script can be run in sbatch mode if needed which requires provision of extra sbatch parameters.
#'
#' @param mutation_file Path to tab-delimited file with columns: tumor_cohort, mutation_class, mutation_file, optional read_depth_file (required)
#' @param ge_dir       Directory to GRanges RDS from parse_genomic_elements() (required)
#' @param result_dir FishHook output directory (required)
#' @param repli_path    Path to replication timing GRanges RDS file (optional)
#' @param dinuc_path    Path to dinucleotide percentage GRanges RDS file (optional)
#' @param trinuc_path   Path to trinucleotide percentage GRanges RDS file (optional)
#' @param id_cap        Maximum number of mutations allowed per sample per region (idcap in master fishHook) 
#' @param id_column     Sample ID column in mutation file (idcol in master fishHook) 
#' @param num_tasks     Number of processes (default 1)
#' @param num_cores     Number of cores (mc.cores in master fishHook)
#' @param tile_size     Tile size
#' @return              Invisibly writes the fishHook reports within [result_dir]
#' @export
run_fishhook <- function(
  mutation_file,
  ge_dir,
  result_dir = NULL,
  repli_path = NULL,
  dinuc_path = NULL,
  trinuc_path = NULL,
  id_cap = 1e9,
  id_column = NULL,
  num_tasks = 1,
  num_cores = 1,
  tile_size = NULL) {

    # check for dependencies
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Please install the 'data.table' package")
    }

    if (!requireNamespace("parallel", quietly = TRUE)) {
        stop("Please install the 'parallel' package")
    }

    # read genomic elements
    if (is.null(ge_dir) || !file.exists(ge_dir)){
        stop("Directory with elements missing")
    } else {
        ge_files <- list.files(ge_dir, pattern="\\.Rds$", full.names=TRUE)
        if (length(ge_files) == 0){
            stop("No RDS files found in directory with elements")
        } else {
            ge_names <- toupper(tools::file_path_sans_ext(basename(ge_files)))
            ge.paths <- setNames(ge_files, ge_names)
        }
    }

    # read mutation file
    mut.df <- data.table::fread(mutation_file)
    if (ncol(mut.df) == 3){
        colnames(mut.df) <- c("tumor", "class", "mut_path")
    } else if (ncol(mut.df) == 4){
        colnames(mut.df) <- c("tumor", "class", "mut_path", "ave_rd_path")
    } else {
        stop("Mutation file must have exactly 3 or 4 columns")
    }

    # check covariate files
    if (!is.null(repli_path) && !file.exists(repli_path)){
        warning("Replication timing file missing")
    }

    if (!is.null(dinuc_path) && !file.exists(dinuc_path)){
        warning("Dinucleotide percentage file missing")
    }

    if (!is.null(trinuc_path) && !file.exists(trinuc_path)){
        warning("Trinucleotide percentage file missing")
    }

    # check if num_cores and id_cap are numeric
    if (!is.numeric(num_cores) || num_cores < 1) {
        stop("num_cores must be a positive integer")
    }
    if (!is.numeric(id_cap) || id_cap <= 0) {
        stop("id_cap must be a positive integer")
    }

    # check if id_column specified
    if (is.null(id_column) || !nzchar(id_column)){
        stop("id_column missing")
    }

    # create fishHook output directory
    if (is.null(result_dir)){
        stop("Please provide output directory to store fishHook results")
    } 
    if (!is.null(result_dir) && !dir.exists(result_dir)){
        dir.create(result_dir, recursive=TRUE, showWarnings=FALSE)
    }

    # define subfunction to run fishHook for each tumor cohort x mutation class x genomic element combination
    job.list <- list()
    for (i in seq_len(nrow(mut.df))){
        tumor <- mut.df$tumor[i] # define tumor cohort
        cls   <- mut.df$class[i] # define mutation class

        # check if mutation file exist
        mutf  <- mut.df$mut_path[i]
        if (!file.exists(mutf)){
            warning("Mutation file at ", mutf, " missing")
            next
        }

        # check if average read depth file exist
        rd <- NULL
        if (ncol(mut.df) == 4) {
            rd_candidate <- mut.df$ave_rd_path[i]
            if (file.exists(rd_candidate)) {
                rd <- rd_candidate
            } else {
                warning("Average read depth file at ", rd_candidate, " missing")
            }
        }

        # loop each element
        for (ge in names(ge.paths)) {
            job.list[[length(job.list) + 1]] <- list(
                tumor = tumor,
                cls = cls,
                mutf = mutf,
                rd  = rd,
                ge = ge,
                ge_path = ge.paths[[ge]])
        }
    }

    run_fishhook_trio <- function(job) {

        message(sprintf("Running FishHook: Tumor = %s | Class = %s | Element = %s", job$tumor, job$cls, job$ge))

        # define subdirectory to store fishHook results
        out.base <- file.path(result_dir, job$tumor, job$cls)
        if (!file.exists(out.base)){
            dir.create(out.base, recursive=TRUE, showWarnings=FALSE)
        }
        out.file <- file.path(out.base, paste0(job$tumor, "_", job$cls, "_", job$ge, "_fishhook_output.txt"))

        # build fishHook command
        args <- c(
        "Rscript",
        file.path("R", "Execute_fishHook.R"),
        "--tiles-file", job$ge_path,
        "--mutation-file", job$mutf,
        "--id-cap", id_cap,
        "--id-column", id_column,
        "--num-cores", num_cores,
        "--replication-time", repli_path,
        "--dinuc", dinuc_path,
        "--trinuc", trinuc_path,
        if (!is.null(job$rd)) c("--read-depth", job$rd) else NULL,
        if (!is.null(tile_size)) c("--tile-size", tile_size) else NULL,
        "--output-file", out.file)

        cmd <- paste(shQuote(args), collapse = " ")
        message("Executing command: ", cmd)
        status <- system(cmd)
        if (status != 0) warning("FishHook failed for: ", out.file)
    }

    # Execute fishHook
    if (num_tasks > 1) {
        parallel::mclapply(job.list, run_fishhook_trio, mc.cores = num_tasks)
    } else {
        lapply(job.list, run_fishhook_trio)
    }

    invisible(NULL)
}


#' @title effector_classifier_byrank
#' @description
#'
#' Takes input fishHook results with recorded excess number of mutations per element and GRanges RDS mutation files. Mutation files must be included 
#' in text file with format: column I: tumor cohort (i.e.tumor type), column II: mutation class (SNV or INDEL), column III: full path to RDS GRanges mutations, 
#' column IV: path to RDS GRanges, average read depth (will be ignored if specified). Function runs ranking of mutations 
#' in each tumor cohort x mutation class by scaled CADD score (CADD_scaled in file), element-specific binary features, significant clustering (MutClust in file),
#' number of mutation-calling tools (NUM_TOOLS in file) and linguistic complexity (LC in file). Element-specific binary features assume having different scores in each GRanges RDS file.
#' Top-ranked mutations matching excess number of mutations in each genomic element
#' as recorded by fishHook and satisfying scaled CADD score >= min_cadd (default 20) and/or positive element-specific binary feature score 
#' are flagged as effector.
#'
#' @param fishhook_dir  Path to the directory containing fishHook output files (required)
#' @param mutation_file Path to tab-delimited file with columns: tumor_cohort, mutation_class, path_to_mutation_file (latest must have columns REF, ALT, CADD_scaled, MutClust, NUM_TOOLS and LC) (required)
#' @param result_dir    Path where each GRanges RDS with effector_flag will be written (required)
#' @param min_cadd      Minimum scaled CADD score for a mutation to be considered (default 20)
#' @param num_tasks     Number of processes (default 1)
#' @param ranking_types Vector of genomic element types to be used for ranking (default ALL)
#' @param coding_types Vector of coding genomic element types to be used for ranking; GRanges RDS must have column LOFTEE (HC/NC) (optional)
#' @param regulatory_types Vector of gene regulatory genomic element types to be used for ranking; GRanges RDS must have columns motif_break_status and motif_gain_status (0 or 1) (optional)
#' @param mirna_binding_types Vector of miRNA-binding genomic element types to be used for ranking; GRanges RDS must have column mirBS_status (0 or 1) (optional)
#' @param mirna_type Label of miRNA element; GRanges RDS must have column rnafold_status (0 or 1) (optional)
#' @param ctcf_type Label of CTCF-binding motifs; GRanges RDS must have columns motif_break_status and motif_gain_status (0 or 1) and motif_name_break (optional)
#' @return Invisibly writes each GRanges mutation file with effector_flag_byrank within [result_dir]
#' @export
effector_classifier_byrank <- function(
  fishhook_dir,
  mutation_file,
  result_dir,
  num_tasks = 1,
  min_cadd = 20,
  ranking_types = "ALL",
  coding_types = NULL,
  regulatory_types = NULL,
  mirna_binding_types = NULL,
  mirna_type = NULL,
  ctcf_type = NULL) {

    # check for dependencies
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Please install the 'data.table' package")
    }
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop("Please install the 'GenomicRanges' package")
    }
    if (!requireNamespace("IRanges", quietly = TRUE)) {
        stop("Please install the 'IRanges' package")
    }

    # define helper functions
    move_to_front <- function(idx, flag) {
        c(idx[flag], idx[!flag])
    }

    binary_positive <- function(idx_vec, seg_type, mcols_df) {

        # check if coding type(s) was/were specified
        if (!is.null(coding_types) && seg_type %in% coding_types) {
            if ("LOFTEE" %in% colnames(mcols_df)) {
                return(mcols_df$LOFTEE[idx_vec] == "HC")
            } else {
                warning("Given following coding types ", coding_types, " but no LOFTEE column found. Skipping ...")
                return(rep(FALSE, length(idx_vec)))
            }
        }

        # check if gene regulatory type(s) was/were specified
        if (!is.null(regulatory_types) && seg_type %in% regulatory_types) {
            if (all(c("motif_gain_status", "motif_break_status") %in% colnames(mcols_df))) {
                return(mcols_df$motif_gain_status[idx_vec] == 1 | mcols_df$motif_break_status[idx_vec] == 1)
            } else {
                warning("Given following gene regulatory types ", regulatory_types, " but no motif_gain_status and/or motif_break_status column found. Skipping ...")
                return(rep(FALSE, length(idx_vec)))
            }
        }

        # check if miRNA-binding type(s) was/were specified
        if (!is.null(mirna_binding_types) && seg_type %in% mirna_binding_types) {
            if ("mirBS_status" %in% colnames(mcols_df)) {
                return(mcols_df$mirBS_status[idx_vec] == 1)
            } else {
                warning("Given following miRNA-binding types ", mirna_binding_types, " but no mirBS_status column found. Skipping ...")                
                return(rep(FALSE, length(idx_vec)))
            }
        }

        # check if miRNA type was specified
        if (!is.null(mirna_type) && seg_type == mirna_type) {
            if ("rnafold_status" %in% colnames(mcols_df)) {
                return(mcols_df$rnafold_status[idx_vec] == 1)
            } else {
                warning("miRNA type specified but no rnafold_status column found. Skipping ...")                
                return(rep(FALSE, length(idx_vec)))
            }
        }

        # CTCF: check motif_gain_status, motif_break_status, motif_name_break
        if (!is.null(ctcf_type) && seg_type == ctcf_type) {
            if (all(c("motif_gain_status", "motif_break_status", "motif_name_break") %in% colnames(mcols_df))) {
                return(mcols_df$motif_gain_status[idx_vec] == 1 | (mcols_df$motif_break_status[idx_vec] == 1 & mcols_df$motif_name_break[idx_vec] == "CTCF"))
            } else {
                warning("CTCF type specified but no motif_gain_status and/or motif_break_status and/or motif_name_break column found. Skipping ...")                
                return(rep(FALSE, length(idx_vec)))
            }
        }

        # If none of the above, or required columns missing, return all FALSE
        return(rep(FALSE, length(idx_vec)))
    }

    # rank candidate indices according to the active criteria
    rank_segment <- function(mut_idx, seg_type, n_excess, mcols_df, use_binary) {
        if (length(mut_idx) == 0 || n_excess <= 0) return(integer())

        required_cols <- c("CADD_scaled", "MutClust", "NUM_TOOLS", "LC")
        missing <- setdiff(required_cols, colnames(mcols_df))
        if (length(missing) > 0) {
            stop("Missing required columns in mutation GRanges: ", paste(missing, collapse = ", "))
        }

        # 1) rank by descending scaled CADD
        ordered_idx <- mut_idx[order(-mcols_df$CADD_scaled[mut_idx], na.last = TRUE)]

        if (!use_binary) {

            # 2) rank by significant clustering
            ordered_idx <- move_to_front(ordered_idx, mcols_df$MutClust[ordered_idx] == 1)

            # 3) rank by descending NUM_TOOLS
            vals_nt <- mcols_df$NUM_TOOLS[ordered_idx]
            for (v in sort(unique(vals_nt), decreasing = TRUE)) {
                ordered_idx <- move_to_front(ordered_idx, vals_nt == v)
            }

            # 4) rank by descending LC
            vals_lc <- mcols_df$LC[ordered_idx]
            for (v in sort(unique(vals_lc), decreasing = TRUE)) {
                ordered_idx <- move_to_front(ordered_idx, vals_lc == v)
            }

            return(head(ordered_idx, min(n_excess, length(ordered_idx))))
        }

        # if use_binary = TRUE, incorporate binary score
        # 2) binary feature
        ordered_idx <- move_to_front(ordered_idx,binary_positive(ordered_idx, seg_type, mcols_df))

        # 3) rank by significant clustering
        ordered_idx <- move_to_front(ordered_idx, mcols_df$MutClust[ordered_idx] == 1)

        # 4) rank by descending NUM_TOOLS
        vals_nt <- mcols_df$NUM_TOOLS[ordered_idx]
        for (v in sort(unique(vals_nt), decreasing = TRUE)) {
            ordered_idx <- move_to_front(ordered_idx, vals_nt == v)
        }

        # 5) rank by descending LC
        vals_lc <- mcols_df$LC[ordered_idx]
        for (v in sort(unique(vals_lc), decreasing = TRUE)) {
            ordered_idx <- move_to_front(ordered_idx, vals_lc == v)
        }

        return(head(ordered_idx, min(n_excess, length(ordered_idx))))
    }

    # read mutation file
    mut.df <- data.table::fread(mutation_file)
    if (ncol(mut.df) < 3) {
        stop("Mutation file must have at least 3 columns")
    } else {
        data.table::setnames(mut.df, 1:3, c("tumor", "class", "mut_path"))
    }

    # create results directory
    if (is.null(result_dir)){
        stop("Output directory to store ranking results missing")
    } 

    if (!is.null(result_dir) && !dir.exists(result_dir)){
        dir.create(result_dir, recursive=TRUE, showWarnings=FALSE)
    }

    # read fishHook output files
    fish_files <- list.files(path = fishhook_dir, pattern = "\\.txt(\\.gz)?$", recursive = TRUE, full.names = TRUE)

    if (length(fish_files) == 0) {
        stop("No FishHook files found in ", fishhook_dir)
    }

    # loop through tumor cohort x mutation class x genomic element and rank mutations
    for (i in seq_len(nrow(mut.df))) {
        tumor <- mut.df$tumor[i] # define tumor cohort
        mut.class   <- mut.df$class[i] # define mutation class
        
        message("Running ranking of elements for tumor: ", tumor, " and mutation class: ", mut.class)

        # detect mode (snv or indel)
        mode_lower <- tolower(mut.class)
        if (!mode_lower %in% c("snv", "indel")) {
            stop("Mutation class must be either 'SNV' or 'INDEL'; got: ", mut.class)
        }

        # check if mutation file exist and if GRanges
        mutf  <- mut.df$mut_path[i]
        if (!file.exists(mutf)){
            warning("Mutation file at ", mutf, " missing")
            next
        }

        mut.gr <- readRDS(mutf)
        if (!inherits(mut.gr, "GRanges")) {
            stop("Object in ", mut.gr, " is not a GRanges")
        }

        # set new metadata column
        mcols(mut.gr)$effector_flag_byrank <- integer(length(mut.gr))

        # check if all mutations are snv or indel
        if (all(c("REF", "ALT") %in% colnames(mcols(mut.gr)))) {
            ref_len     <- nchar(as.character(mcols(mut.gr)$REF))
            alt_max_len <- vapply(
                mcols(mut.gr)$ALT,
                function(x) max(nchar(as.character(x))),
                integer(1L))
            is_snv_flag <- ref_len == 1 & alt_max_len == 1
            wrong <- if (mode_lower == "snv") which(!is_snv_flag) else which(is_snv_flag)
            if (length(wrong) > 0) {
                stop("Mutation GRanges contains ", length(wrong), if (mode_lower == "snv") " INDEL(s)" else " SNV(s)", " but mutation_class = '", mut.class, "'.")
            }
        } else {
            stop("REF and/or ALT columns not found")
        }

        # select corresponding fishHook output for tumor cohort x mutation class
        pattern_tc <- paste0("/", tumor, "/", mut.class, "/")
        fh_subset <- fish_files[grepl(pattern_tc, fish_files, fixed = TRUE)]
        if (!length(fh_subset)) next

        # select all ranking types
        if (identical(ranking_types, "ALL")) {
            all_types <- character()
            for (fh in fh_subset) {
                dt_tmp   <- data.table::fread(fh, select = "type")
                all_types <- unique(c(all_types, dt_tmp$type))
            }
            allowed_types_eff <- all_types
        } else {
            allowed_types_eff <- ranking_types
        }

        # loop through each fishHook output
        for (fh in fh_subset) {
            dt <- data.table::fread(fh)
            required_cols <- c("seqnames", "start", "end", "type", "obs_minus_exp")
            if (!all(required_cols %in% names(dt))) {
                stop("Missing columns in FishHook file: ", fh, "\n Required: ", paste(required_cols, collapse = ", "))
            }

            # keep only rows with obs_minus_exp > 0 and allowed types
            dt <- dt[obs_minus_exp > 0 & type %in% allowed_types_eff]

            if (nrow(dt) == 0) next
            # convert fishHook output to GRanges and find overlapping mutations 
            seg_gr <- GenomicRanges::GRanges(
                seqnames = as.character(dt$seqnames),
                ranges   = IRanges::IRanges(as.integer(dt$start), as.integer(dt$end)),
                strand   = "*",
                type     = as.character(dt$type),
                n_excess = as.integer(round(dt$obs_minus_exp)))

            ov <- GenomicRanges::findOverlaps(seg_gr, mut.gr, ignore.strand = TRUE)
            if (!length(ov)) next
            for (sid in unique(queryHits(ov))) {
                seg_type <- mcols(seg_gr)$type[sid] # element type
                n_excess <- mcols(seg_gr)$n_excess[sid] # excess mutations
                seg_mut_idx <- subjectHits(ov)[
                    queryHits(ov) == sid]

                # exclude mutations already flagged
                avail_idx <- seg_mut_idx[mcols(mut.gr)$effector_flag_byrank[seg_mut_idx] == 0L]

                # CADD/binary score filter
                is_binary_type <- (
                    (!is.null(coding_types) && seg_type %in% coding_types) ||
                    (!is.null(regulatory_types) && seg_type %in% regulatory_types) ||
                    (!is.null(mirna_binding_types) && seg_type %in% mirna_binding_types) ||
                    (!is.null(mirna_type) && seg_type == mirna_type) ||
                    (!is.null(ctcf_type) && seg_type == ctcf_type))

                if (is_binary_type) {
                    keep <- (mcols(mut.gr)$CADD_scaled[avail_idx] >= min_cadd) | binary_positive(avail_idx, seg_type, mcols(mut.gr))
                } else {
                    keep <- (mcols(mut.gr)$CADD_scaled[avail_idx] >= min_cadd)
                }
                
                cand_idx <- seg_mut_idx[keep]
                if (length(cand_idx) == 0) next

                # subtract flagged mutations from excess
                already_flagged_count <- length(seg_mut_idx) - length(avail_idx)
                n_needed <- n_excess - already_flagged_count
                if (n_needed <= 0) next

                # rank mutations in each element
                chosen <- rank_segment(
                    mut_idx = cand_idx,
                    seg_type = seg_type,
                    n_excess = n_needed,
                    mcols_df = mcols(mut.gr),
                    use_binary = is_binary_type)
                if (length(chosen)) {
                    mcols(mut.gr)$effector_flag_byrank[chosen] <- 1L
                }
            }
        }

        # save as GRanges RDS with effector_flag_byrank column
        out_file <- file.path(result_dir, paste0(tumor, "_", mut.class, "_effector_flag_byrank.Rds"))
        saveRDS(mut.gr, out_file)
    }

    # if (num_tasks > 1) {
    #     # use mclapply to fork [num_tasks] workers, each handling distinct GRanges RDS mutation file
    #     parallel::mclapply(seq_len(nrow(mut.df)), FUN = run_ranking_trio, mc.cores  = num_tasks, mc.silent = TRUE)
    # } else {
    #     # serial execution if [num_tasks] = 1
    #     for (i in seq_len(nrow(mut.df))) {
    #         run_ranking_trio(i)
    #     }
    # }
    
    invisible(NULL)
}


#' @title effector_classifier_byrule
#' @description
#'
#' Takes input GRanges RDS mutation files from ranking approach and parsed genomic element files and assigns an effector flag on overlapping mutations per element type
#' based on the rules in SNV or INDEL mode. 
#' In SNV mode: I) if element type = coding type, call an effector if 1) nonsense and/or start loss mutation with LOFTEE = HC,
#' 2) missense mutation with CADD_scaled >= 20. II) if element type = regulatory type, call an effector if motif_break_status and/or motif_gain_status = 1 & NUM_TOOLS >= 4
#' & LC >= 15. III) if element type = miRNA-binding type, call an effector if mirBS_status = 1 & NUM_TOOLS >= 4 & LC >= 15. IV) if element type = miRNA type,
#' call an effector if rnafold_status = 1. V) if element type = intron type, call an effector if spliceSiteMod_status = 1 & CADD_scaled >= 20. VI) if element type = CTCF type,
#' call an effector if motif_gain_status = 1 & NUM_TOOLS >= 4 & LC >= 15, 2) motif_break_status = 1 & motif_name_break = CTCF & NUM_TOOLS >= 4 & LC >= 15.
#' In INDEL mode: I) if element type = coding type, call an effector if 1) in-frame mutation with CADD_scaled >= 20, 2) frameshift mutation with LOFTEE = HC.
#' II) if element type = regulatory type, call an effector if motif_break_status and/or motif_gain_status = 1 & NUM_TOOLS >= 3 & LC >= 15. 
#' III) if element type = miRNA-binding type, call an effector if mirBS_status = 1 & NUM_TOOLS >= 3 & LC >= 15. IV) if element type = miRNA type,
#' call an effector if rnafold_status = 1. V) if element type = intron type, call an effector if spliceSiteMod_status = 1 & CADD_scaled >= 20. VI) if element type = CTCF type,
#' call an effector if motif_gain_status = 1 & NUM_TOOLS >= 3 & LC >= 15, 2) motif_break_status = 1 & motif_name_break = CTCF & NUM_TOOLS >= 3 & LC >= 15.
#' In both SNV and INDEL modes, for each element type if CADD_scaled exceeds 95th percentile of all CADD_scaled scores within that element type & CADD_scaled >= 35, 
#' call an effector.
#'
#' @param input_rds Input GRanges RDS mutation file (required)
#' @param ge_dir  Directory to GRanges RDS from parse_genomic_elements() (required)
#' @param result_dir Path where each GRanges RDS with effector_flag_byrule will be written (default dirname(input_rds)) 
#' @param mode Mutation class to label as effector or non-effector 
#' @param coding_types Vector of coding genomic element types to be used for ranking; GRanges RDS must have column LOFTEE (HC/NC) and CONSEQUENCE (e.g. missense) (optional)
#' @param regulatory_types Vector of gene regulatory genomic element types to be used for ranking; GRanges RDS must have columns motif_break_status and motif_gain_status (0 or 1) (optional)
#' @param mirna_binding_types Vector of miRNA-binding genomic element types to be used for ranking; GRanges RDS must have column mirBS_status (0 or 1) (optional)
#' @param mirna_type Label of miRNA element; GRanges RDS must have column rnafold_status (0 or 1) (optional)
#' @param ctcf_type Label of CTCF-binding motifs; GRanges RDS must have columns motif_break_status and motif_gain_status (0 or 1) and motif_name_break (optional)
#' @param intron_types Vector of intron elements; GRanges RDS must have column spliceSiteMod_status (0 or 1) (optional)
#'
#' @return Invisibly writes each GRanges mutation file with effector_flag_byrule within [result_dir]
#' @export
effector_classifier_byrule <- function(
    input_rds,
    ge_dir,
    result_dir = dirname(input_rds),
    mode = c("SNV", "INDEL"),
    coding_types = NULL,
    regulatory_types = NULL,
    mirna_binding_types = NULL,
    mirna_type = NULL,
    ctcf_type = NULL,
    intron_types = NULL) {

    mode <- match.arg(mode)

    # read genomic elements
    if (is.null(ge_dir) || !file.exists(ge_dir)){
        stop("Directory with elements missing")
    } else {
        ge_files <- list.files(ge_dir, pattern="\\.Rds$", full.names=TRUE)
        if (length(ge_files) == 0){
            stop("No RDS files found in directory with elements")
        } else {
            ge_names <- toupper(tools::file_path_sans_ext(basename(ge_files)))
            ge.paths <- setNames(ge_files, ge_names)
        }
    }
    types <- tolower(names(ge.paths))

    # read input GRanges RDS mutation file
    mut.gr <- readRDS(input_rds)

    # extract auxiliary columns 
    mcols_df <- as.data.frame(mcols(mut.gr))

    # check if required columns exist
    required_cols <- "CADD_scaled"
    if (!is.null(coding_types)){
        coding_types <- tolower(coding_types)
        missing <- setdiff(coding_types, types)
        if (length(missing) == 0L){
            required_cols <- c(required_cols, "LOFTEE", "CONSEQUENCE")
        } else if (length(missing) > 0 && length(missing) < length(coding_types)) {
            required_cols <- c(required_cols, "LOFTEE", "CONSEQUENCE")
            warning("Following types not found in elements file: ", paste(missing, collapse = ", "), " Omitting...")
        } else {
            warning("All specified coding types not found in elements file.")
            coding_types <- NULL
        }
    }

    if (!is.null(regulatory_types)) {
        regulatory_types <- tolower(regulatory_types)
        missing <- setdiff(regulatory_types, types)
        if (length(missing) == 0L){
            required_cols <- c(required_cols, "motif_break_status", "motif_gain_status", "NUM_TOOLS", "LC")
        } else if (length(missing) > 0 && length(missing) < length(regulatory_types)) {
            required_cols <- c(required_cols, "motif_break_status", "motif_gain_status", "NUM_TOOLS", "LC")
            warning("Following types not found in elements file: ", paste(missing, collapse = ", "), " Omitting...")
        } else {
            warning("All specified regulatory types not found in elements file.")
            regulatory_types <- NULL
        }
    }


    if (!is.null(mirna_binding_types)) {
        mirna_binding_types <- tolower(mirna_binding_types)
        missing <- setdiff(mirna_binding_types, types)
        if (length(missing) == 0L){
            required_cols <- c(required_cols, "mirBS_status", "NUM_TOOLS", "LC")
        } else if (length(missing) > 0 && length(missing) < length(mirna_binding_types)) {
            required_cols <- c(required_cols, "mirBS_status", "NUM_TOOLS", "LC")
            warning("Following types not found in elements file: ", paste(missing, collapse = ", "), " Omitting...")
        } else {
            warning("All specified miRNA-binding types not found in elements file.")
            mirna_binding_types <- NULL
        }
    }


    if (!is.null(mirna_type)) {
        mirna_type <- tolower(mirna_type)
        missing <- setdiff(mirna_type, types)
        if (length(missing) == 0L){
            required_cols <- c(required_cols, "rnafold_status")
        } else {
            warning("Specified miRNA not found in elements file.")
            mirna_type <- NULL
        }
    }


    if (!is.null(intron_types)) {
        intron_types <- tolower(intron_types)
        missing <- setdiff(intron_types, types)
        if (length(missing) == 0L){
            required_cols <- c(required_cols, "spliceSiteMod_status")
        } else if (length(missing) > 0 && length(missing) < length(intron_types)) {
            required_cols <- c(required_cols, "spliceSiteMod_status")
            warning("Following types not found in elements file: ", paste(missing, collapse = ", "), " Omitting...")
        } else {
            warning("All specified intron types not found in elements file.")
            intron_types <- NULL
        }
    }


    if (!is.null(ctcf_type)) {
        ctcf_type <- tolower(ctcf_type)
        missing <- setdiff(ctcf_type, types)
        if (length(missing) == 0L){
            required_cols <- c(required_cols, "motif_break_status", "motif_gain_status", "NUM_TOOLS", "LC", "motif_name_break")
        } else {
            warning("Specified CTCF type not found in elements file.")
            ctcf_type <- NULL
        }
    }
    

    missing_cols <- setdiff(unique(required_cols), names(mcols(mut.gr)))
    if (length(missing_cols) > 0L) {
        stop("Missing required metadata columns in mutation GRanges: ", paste(missing_cols, collapse = ", "))
    }

    # define effector_flag vector
    nmut <- length(mut.gr)
    effector_flag <- rep(FALSE, nmut)

    # assign required criteria for each mutation within element
    process_element <- function(elem.gr, type, nt.threshold) {

        # get overlapping mutations with element type
        hits <- findOverlaps(mut.gr, elem.gr, ignore.strand = TRUE)
        if (length(hits) == 0L) return(invisible(NULL))
        hits.q <- queryHits(hits)
        unique_idx <- unique(hits.q)

        # extract 95th percentile CADD_scaled for these mutations
        cadd_vals <- mcols_df$CADD_scaled[unique_idx]
        q95 <- if (all(is.na(cadd_vals))) NA else as.numeric(quantile(cadd_vals, 0.95, na.rm = TRUE))

        # define required columns
        consequence <- mcols_df$CONSEQUENCE[unique_idx]
        loftee <- mcols_df$LOFTEE[unique_idx]
        motif_break_status <- mcols_df$motif_break_status[unique_idx]
        motif_gain_status <- mcols_df$motif_gain_status[unique_idx]
        num_tools <- mcols_df$NUM_TOOLS[unique_idx]
        lc <- mcols_df$LC[unique_idx]
        mirbs_status <- mcols_df$mirBS_status[unique_idx]
        rnafold_status <- mcols_df$rnafold_status[unique_idx]
        splice_status <- mcols_df$spliceSiteMod_status[unique_idx]
        motif_name_break <- mcols_df$motif_name_break[unique_idx]
        cadd <- cadd_vals

        eff.flag <- logical(length(unique_idx))

        # 1) Rules for coding type elements
        if (type %in% coding_types) {
            if (mode == "SNV") {
                eff.flag <- (grepl("\\bnonsense\\b|\\bstart_loss\\b", consequence) & loftee == "HC") |
                            (grepl("\\bmissense\\b", consequence) & !is.na(cadd) & cadd >= 20)
            } else {
                is_inframe <- grepl("inframe", consequence, ignore.case = TRUE)
                is_frameshift <- grepl("frameshift", consequence, ignore.case = TRUE)
                eff.flag <- (is_inframe & !is.na(cadd) & cadd >= 20) |
                            (is_frameshift & loftee == "HC")
            }
        }

        # 2) Rules for regulatory elements
        if (type %in% regulatory_types) {
            eff.flag <- eff.flag |
                (((motif_break_status == 1L | motif_gain_status == 1L) &
                  !is.na(num_tools) & num_tools >= nt.threshold &
                  !is.na(lc) & lc >= 15))
        }

        # 3) Rules for miRNA-binding elements
        if (type %in% mirna_binding_types) {
            eff.flag <- eff.flag |
                (mirbs_status == 1L & !is.na(num_tools) & num_tools >= nt.threshold &
                 !is.na(lc) & lc >= 15)
        }

        # 4) Rules for miRNA element
        if (!is.null(mirna_type) && type == mirna_type) {
            eff.flag <- eff.flag | (rnafold_status == 1L)
        }

        # 5) Rules for intron elements
        if (type %in% intron_types) {
            eff.flag <- eff.flag |
                (!is.na(splice_status) & splice_status == 1L &
                 !is.na(cadd) & cadd >= 20)
        }

        # 6) Rules for CTCF-binding elements
        if (type %in% ctcf_type) {
            ctcf_gain <- (motif_gain_status == 1L & !is.na(num_tools) & num_tools >= nt.threshold & !is.na(lc) & lc >= 15)
            ctcf_break <- (motif_break_status == 1L & motif_name_break == "CTCF" & !is.na(num_tools) & num_tools >= nt.threshold & !is.na(lc) & lc >= 15)
            eff.flag <- eff.flag | ctcf_gain | ctcf_break
        }

        # 7) Universal CADD rule for this element
        eff.flag <- eff.flag | (!is.na(cadd) & cadd >= 35 & !is.na(q95) & cadd >= q95)
        effector_flag[unique_idx[eff.flag]] <<- TRUE

        invisible(NULL)
    }

    # process each provided element category
    i <- 0
    for (type.upper in names(ge.paths)) {
        elem.path <- ge.paths[[type.upper]]
        type <- tolower(type.upper)
        elem.gr <- readRDS(elem.path)

        # set threshold for NUM_TOOLS based on mode
        thr <- if (mode == "SNV") 4 else 3

        i = i + 1
        message(sprintf("Processing [%d/%d]: %s", i, length(ge.paths), type.upper))

        # process each element
        if (type %in% coding_types) {
            process_element(elem.gr, type, thr)
        } else if (type %in% regulatory_types) {
            process_element(elem.gr, type, thr)
        } else if (type %in% mirna_binding_types) {
            process_element(elem.gr, type, thr)
        } else if (type == mirna_type) {
            process_element(elem.gr, type, NA)
        } else if (type %in% intron_types) {
            process_element(elem.gr, type, NA)
        } else if (type %in% ctcf_type) {
            process_element(elem.gr, type, thr)
        } else {
            process_element(elem.gr, type, NA)
        }
    }

    message("Finished classification. Total mutations flagged: ", sum(effector_flag))

    # attach effector_flag_byrule to mutation GRanges
    mcols(mut.gr)$effector_flag_byrule <- as.integer(effector_flag)

    # save output as RDS file
    if (!dir.exists(result_dir)) {
        dir.create(result_dir, recursive = TRUE)
    }

    outpath <- file.path(result_dir, paste0(sub("\\.Rds$", "", basename(input_rds)), "_effectorflag_byrule.Rds"))
    saveRDS(mut.gr, outpath)

    invisible(NULL)
}
