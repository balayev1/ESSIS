############################################################################
## essis.R
## Effector Somatic SNV and INDEL Seeker (ESSIS)
############################################################################


#' @title essis
#' @description
#'
#' Split into six steps: 1) parses provided genomic elements into list of GRanges 
#' 2) runs fishHook using each GRanges element set as hypotheses and each mutation file representative of each tumor cohort x mutation class as events.
#' Mutation classes are single-nucleotide variant (SNV) and insertion/deletion (INDEL). Each mutation file MUST contain columns "CADD_scaled", "NUM_TOOLS", "LC",
#' "CONSEQUENCE", "SAMPLEID", "REF", "ALT" and other optional columns depending on the specified element types given as input. Covariates such as replication timing, dinucleotide percentages, 
#' trinucleotide percentages, average read depth and local mutation density can be included within fishHook model. 
#' Excess number of mutations is reported per each genomic element per each tumor cohort x mutation class.
#' 3) runs ranking of mutations in each tumor cohort x mutation class by scaled CADD score, element-specific binary features, significant clustering,
#' number of mutation-calling tools and linguistic complexity. Top-ranked mutations matching excess number of mutations in each genomic element
#' as recorded by fishHook are flagged as effector. 
#' 4) runs fishHook using GRanges element set as hypotheses and mutations merged across all tumor cohorts by mutation class as events which also excluded 
#' flagged effector mutations from previous step. 
#' 5) runs ranking of mutations in each merged mutation file by same criteria as in first ranking step and flags top-ranked mutations matching excess number 
#' of mutations in each genomic element from fishHook as effectors.
#' 6) runs rule-based approach on all mutations failing the effector flag from previous steps. Please see following page for the rules ().
#'
#' @param ge_file   Full path to GRanges RDS with genomic elements (required)
#' @param gtf_path  Full path to a GTF annotation to acquire intergenic regions (optional)
#' @param output_dir Directory to store all pipeline results (required)
#' @param mut_file Path to tab-delimited file with columns: tumor_cohort, mutation_class, mutation_file, optional read_depth_file (required)
#' @param repli Path to replication timing GRanges RDS file (optional)
#' @param dinuc Path to dinucleotide percentage GRanges RDS file (optional)
#' @param trinuc Path to trinucleotide percentage GRanges RDS file (optional)
#' @param id_cap Maximum number of mutations allowed per sample per region (idcap in master fishHook) 
#' @param id_col Sample ID column in mutation file (idcol in master fishHook) 
#' @param tasks Number of processes (default 1)
#' @param cores Number of cores (mc.cores in master fishHook)
#' @param tile_size Tile size
#' @param min_cadd Minimum scaled CADD score for a mutation to be considered (default 20)
#' @param ranking_types Vector of genomic element types to be used for ranking (default ALL)
#' @param coding_types Vector of coding genomic element types to be used for ranking; GRanges RDS must have column LOFTEE (HC/NC) (optional)
#' @param regulatory_types Vector of gene regulatory genomic element types to be used for ranking; GRanges RDS must have columns motif_break_status and motif_gain_status (0 or 1) (optional)
#' @param mirna_binding_types Vector of miRNA-binding genomic element types to be used for ranking; GRanges RDS must have column mirBS_status (0 or 1) (optional)
#' @param mirna_type Label of miRNA element; GRanges RDS must have column rnafold_status (0 or 1) (optional)
#' @param ctcf_type Label of CTCF-binding motifs; GRanges RDS must have columns motif_break_status and motif_gain_status (0 or 1) and motif_name_break (optional)
#' @param intron_types Vector of intron elements; GRanges RDS must have column spliceSiteMod_status (0 or 1) (optional)
#'
#' @return  GRanges RDS with effector flag and all other included variables from input GRanges RDS files
#' @export

# check for dependencies
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Please install the 'GenomicRanges' package")
}

if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Please install the 'rtracklayer' package")
}

if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("Please install the 'optparse' package")
}

options(warn = -1)
suppressPackageStartupMessages({library(GenomicRanges)})
suppressPackageStartupMessages({library(rtracklayer)})
suppressPackageStartupMessages({library(optparse)})

# source the file with functions
suppressWarnings(suppressMessages(suppressPackageStartupMessages(source(file.path("R", "functions.R")))))

# define parsed arguments
option_list <- list(
    make_option("--ge-file", type = "character", help = "Path to input RDS file with genomic elements (required)"),
    make_option("--gtf-path", type = "character", default=NULL, help = "Path to GTF file"),
    make_option("--output-dir", type = "character", default = NULL, help = "Output directory to store all results (required)"),
    make_option("--mut-file", type="character", help="Tab-delimited file: tumor_cohort, mutation_class, mutation_file_path (in RDS format), optional ave_rd_path (in RDS format) (required)"),
    make_option("--repli", type="character", default=NULL, help="Replication timing RDS (optional)"),
    make_option("--dinuc", type="character", default=NULL, help="Dinucleotide % RDS (optional)"),
    make_option("--trinuc", type="character", default=NULL, help="Trinucleotide % RDS (optional)"),
    make_option("--id-cap", type="numeric", default=1e9, help="idcap in master fishHook"),
    make_option("--id-col", type="character", default=NULL, help="idcol in master fishHook"),
    make_option("--tasks", type="integer", default=1, help="Number of processes (default 1)"),
    make_option("--cores", type="integer", default=1, help="mc.cores in master fishHook"),
    make_option("--tile-size", type="integer", default=NULL, help="split size of genome tiles"),
    make_option("--min-cadd", type = "numeric", default = 20, help = "Minimum scaled CADD score for consideration [default %default]"),
    make_option("--ranking-types", type = "character", default = 'ALL', help = "Comma-separated list of element types to rank (e.g. CDS,3UTR,superenhancer)"),
    make_option("--coding-types", type = "character", default = NULL, help = "Comma-separated list of coding element types (e.g. CDS); requires column LOFTEE in mutation RDS file with factors HC=high-confidence and NC=no-confidence. Please check LOFTEE github page for more."),
    make_option("--regulatory-types", type = "character", default = NULL, help = "Comma-separated list of regulatory element types (e.g. promoters,superenhancers); requires columns motif_gain_status and motif_break_status in mutation RDS file with factors 0 and 1."),
    make_option("--mirna_binding-types", type = "character", default = NULL, help = "Comma-separated list of miRNA-binding element types (e.g. 3UTR); requires column mirBS_status in mutation RDS file with factors 0 and 1."),
    make_option("--mirna-type", type = "character", default = NULL, help = "miRNA label"),
    make_option("--ctcf-type", type = "character", default = NULL, help = "CTCF-binding motif label"),
    make_option("--intron-types", type = "character", default = NULL, help = "Comma-separated list of intronic element types (e.g. intron); requires column spliceSiteMod_status in mutation RDS file with factors 0 and 1."))

opt <- parse_args(OptionParser(option_list = option_list))

# run the pipeline
# Step I: 
message("Parsing genomic elements ...\n")

ge.list <- parse_genomic_elements(
    ge_file = opt$`ge-file`,
    gtf_path = opt$`gtf-path`,
    output_dir = file.path(opt$`output-dir`, "parsed_elements"))

message("Parsing is done.\n")

# Step II:
message("Running fishHook on each tumor cohort x mutation class x genomic element...\n")

# path to output directory from parse_genomic_elements()
ge_dir <- file.path(opt$`output-dir`, "parsed_elements")

run_fishhook(
    mutation_file = opt$`mut-file`,
    ge_dir = ge_dir,
    result_dir = file.path(opt$`output-dir`, "fishHook", "stepI"),
    repli_path = opt$repli,
    dinuc_path = opt$dinuc,
    trinuc_path = opt$trinuc,
    id_cap = opt$`id-cap`,
    id_column = opt$`id-col`,
    num_tasks = opt$tasks,
    num_cores = opt$cores,
    tile_size = opt$`tile-size`)

message("FishHook run is done.\n")

# Step III:
message("Ranking mutations on each tumor cohort x mutation class x genomic element...\n")

effector_classifier_byrank(
    fishhook_dir = file.path(opt$`output-dir`, "fishHook", "stepI"),
    mutation_file = opt$`mut-file`,
    result_dir = file.path(opt$`output-dir`, "ranking", "stepI"),
    num_tasks = opt$tasks,
    min_cadd = opt$`min-cadd`,
    ranking_types = if (identical(opt$`ranking-types`, "ALL")) "ALL" else strsplit(opt$`ranking-types`, ",")[[1]],
    coding_types = if (is.null(opt$`coding-types`)) NULL else strsplit(opt$`coding-types`, ",")[[1]],
    regulatory_types = if (is.null(opt$`regulatory-types`)) NULL else strsplit(opt$`regulatory-types`, ",")[[1]],
    mirna_binding_types = if (is.null(opt$`mirna_binding-types`)) NULL else strsplit(opt$`mirna_binding-types`, ",")[[1]],
    mirna_type = if (nzchar(opt$`mirna-type`)) opt$`mirna-type` else NULL,
    ctcf_type = if (nzchar(opt$`ctcf-type`) ) opt$`ctcf-type` else NULL)

message("Ranking is done.\n")


# Step IV:
message("Running fishHook on all tumor cohorts x mutation class after excluding effector-flagged mutations...\n")

# read SNV and INDEL ranking RDS outputs
intermediate.dir <- file.path(opt$`output-dir`, "intermediate_files")
if (!dir.exists(intermediate.dir)){
    dir.create(intermediate.dir, recursive = TRUE, showWarnings = TRUE)
}

all_rds <- list.files(path = file.path(opt$`output-dir`, "ranking", "stepI"), 
    pattern = "\\_effector_flag_byrank.Rds$", full.names = TRUE)
snv_files   <- all_rds[grepl("_SNV_",   basename(all_rds), fixed = TRUE)]
indel_files <- all_rds[grepl("_INDEL_", basename(all_rds), fixed = TRUE)]

# merge SNV and INDEL files separately after excluding effector-flagged mutations
tumor <- c(); class <- c(); mut_path <- c()
if (length(snv_files)) {
    snv_list <- lapply(snv_files, function(f) {
        gr <- readRDS(f)
        gr[mcols(gr)$effector_flag_byrank == 0]})

    merged_snv <- do.call(c, snv_list)

    saveRDS(merged_snv, file = file.path(intermediate.dir, "merged_SNV_non_effector_from_stepI.Rds"))

    tumor <- append(tumor, 'ALL')
    class <- append(class, 'SNV')
    mut_path <- append(mut_path, file.path(intermediate.dir, "merged_SNV_non_effector_from_stepI.Rds"))
}

if (length(indel_files)) {
    indel_list <- lapply(indel_files, function(f) {
        gr <- readRDS(f)
        gr[mcols(gr)$effector_flag_byrank == 0]})

    merged_indel <- do.call(c, indel_list)

    saveRDS(merged_indel,file = file.path(intermediate.dir, "merged_INDEL_non_effector_from_stepI.Rds"))

    tumor <- append(tumor, 'ALL')
    class <- append(class, 'INDEL')
    mut_path <- append(mut_path, file.path(intermediate.dir, "merged_INDEL_non_effector_from_stepI.Rds"))
}

temp.df <- data.frame(tumor = tumor, class = class, mut_path = mut_path)
if (nrow(temp.df) > 0){
    write.table(temp.df, file = file.path(intermediate.dir, "all_muts.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

    run_fishhook(
        mutation_file = file.path(intermediate.dir, "all_muts.txt"),
        ge_dir = ge_dir,
        result_dir = file.path(opt$`output-dir`, "fishHook", "stepII"),
        repli_path = opt$repli,
        dinuc_path = opt$dinuc,
        trinuc_path = opt$trinuc,
        id_cap = opt$`id-cap`,
        id_column = opt$`id-col`,
        num_tasks = opt$tasks,
        num_cores = opt$cores,
        tile_size = opt$`tile-size`)

    message("FishHook run is done.\n")
} else {
    stop("FishHook run aborted. No input files found ...\n")
}

# subset effectors as well
if (length(snv_files)) {
    snv_list <- lapply(snv_files, function(f) {
    gr <- readRDS(f)
    gr[mcols(gr)$effector_flag_byrank == 1]})

    merged_snv <- do.call(c, snv_list)

    saveRDS(merged_snv, file = file.path(intermediate.dir, "merged_SNV_effector.Rds"))

}

if (length(indel_files)) {
    indel_list <- lapply(indel_files, function(f) {
        gr <- readRDS(f)
        gr[mcols(gr)$effector_flag_byrank == 1]})

    merged_indel <- do.call(c, indel_list)

    saveRDS(merged_indel,file = file.path(intermediate.dir, "merged_INDEL_effector.Rds"))
}


# Step V:
message("Ranking mutations on all tumor cohorts x mutation class\n")

effector_classifier_byrank(
    fishhook_dir = file.path(opt$`output-dir`, "fishHook", "stepII"),
    mutation_file = file.path(intermediate.dir, "all_muts.txt"),
    result_dir = file.path(opt$`output-dir`, "ranking", "stepII"),
    num_tasks = opt$tasks,
    min_cadd = opt$`min-cadd`,
    ranking_types = if (identical(opt$`ranking-types`, "ALL")) "ALL" else strsplit(opt$`ranking-types`, ",")[[1]],
    coding_types = if (is.null(opt$`coding-types`)) NULL else strsplit(opt$`coding-types`, ",")[[1]],
    regulatory_types = if (is.null(opt$`regulatory-types`)) NULL else strsplit(opt$`regulatory-types`, ",")[[1]],
    mirna_binding_types = if (is.null(opt$`mirna_binding-types`)) NULL else strsplit(opt$`mirna_binding-types`, ",")[[1]],
    mirna_type = if (nzchar(opt$`mirna-type`)) opt$`mirna-type` else NULL,
    ctcf_type = if (nzchar(opt$`ctcf-type`) ) opt$`ctcf-type` else NULL)

message("Ranking is done.\n")

# Step VI:
message("Rule-based annotation of mutations across all tumor cohorts...\n")

# combine all effector mutations from previous ranking at step III and step V
final.results.dir <- file.path(opt$`output-dir`, "final_results")
if (!dir.exists(final.results.dir)){
    dir.create(final.results.dir, recursive = TRUE, showWarnings = TRUE)
}

snv.eff.path <- file.path(intermediate.dir, "merged_SNV_effector.Rds")
snv.final.rank.path <- file.path(opt$`output-dir`, "ranking", "stepII", "ALL_SNV_effector_flag_byrank.Rds")
if (file.exists(snv.eff.path) && file.exists(snv.final.rank.path)){
    snv.eff <- readRDS(snv.eff.path)
    snv.final.rank <- readRDS(snv.final.rank.path)
    saveRDS(do.call(c, list(snv.eff, snv.final.rank)), file = file.path(intermediate.dir, "postrank_SNV.Rds"))

    effector_classifier_byrule(
        input_rds = file.path(intermediate.dir, "postrank_SNV.Rds"),
        ge_dir = ge_dir,
        result_dir = final.results.dir,
        mode = "SNV",
        coding_types = if (is.null(opt$`coding-types`)) NULL else strsplit(opt$`coding-types`, ",")[[1]],
        regulatory_types = if (is.null(opt$`regulatory-types`)) NULL else strsplit(opt$`regulatory-types`, ",")[[1]],
        mirna_binding_types = if (is.null(opt$`mirna_binding-types`)) NULL else strsplit(opt$`mirna_binding-types`, ",")[[1]],
        mirna_type = if (nzchar(opt$`mirna-type`)) opt$`mirna-type` else NULL,
        ctcf_type = if (nzchar(opt$`ctcf-type`) ) opt$`ctcf-type` else NULL,
        intron_types = if (is.null(opt$`intron-types`)) NULL else strsplit(opt$`intron-types`, ",")[[1]])
}

indel.eff.path <- file.path(intermediate.dir, "merged_INDEL_effector.Rds")
indel.final.rank.path <- file.path(opt$`output-dir`, "ranking", "stepII", "ALL_INDEL_effector_flag_byrank.Rds")
if (file.exists(indel.eff.path) && file.exists(indel.final.rank.path)){
    indel.eff <- readRDS(indel.eff.path)
    indel.final.rank <- readRDS(indel.final.rank.path)
    saveRDS(do.call(c, list(indel.eff, indel.final.rank)), file = file.path(intermediate.dir, "postrank_INDEL.Rds"))

    effector_classifier_byrule(
        input_rds = file.path(intermediate.dir, "postrank_INDEL.Rds"),
        ge_dir = ge_dir,
        result_dir = final.results.dir,
        mode = "INDEL",
        coding_types = if (is.null(opt$`coding-types`)) NULL else strsplit(opt$`coding-types`, ",")[[1]],
        regulatory_types = if (is.null(opt$`regulatory-types`)) NULL else strsplit(opt$`regulatory-types`, ",")[[1]],
        mirna_binding_types = if (is.null(opt$`mirna_binding-types`)) NULL else strsplit(opt$`mirna_binding-types`, ",")[[1]],
        mirna_type = if (nzchar(opt$`mirna-type`)) opt$`mirna-type` else NULL,
        ctcf_type = if (nzchar(opt$`ctcf-type`) ) opt$`ctcf-type` else NULL,
        intron_types = if (is.null(opt$`intron-types`)) NULL else strsplit(opt$`intron-types`, ",")[[1]])
}

message("Rule-based annotation of mutations is done.\n")


# Step VII:
message("Merging effector flags from rank- and rule-based classification...\n")

snv.muts.gr <- readRDS(file.path(final.results.dir, "postrank_SNV_effectorflag_byrule.Rds"))
rule_flag <- if ("effector_flag_byrule" %in% names(mcols(snv.muts.gr))) mcols(snv.muts.gr)$effector_flag_byrule else 0L
rank_flag <- if ("effector_flag_byrank" %in% names(mcols(snv.muts.gr))) mcols(snv.muts.gr)$effector_flag_byrank else 0L
mcols(snv.muts.gr)$effector_flag <- as.integer((rank_flag == 1L) | (rule_flag == 1L))
saveRDS(snv.muts.gr, file.path(final.results.dir, "postrank_SNV_effectorflagged.Rds"))

indel.muts.gr <- readRDS(file.path(final.results.dir, "postrank_INDEL_effectorflag_byrule.Rds"))
rule_flag <- if ("effector_flag_byrule" %in% names(mcols(indel.muts.gr))) mcols(indel.muts.gr)$effector_flag_byrule else 0L
rank_flag <- if ("effector_flag_byrank" %in% names(mcols(indel.muts.gr))) mcols(indel.muts.gr)$effector_flag_byrank else 0L
mcols(indel.muts.gr)$effector_flag <- as.integer((rank_flag == 1L) | (rule_flag == 1L))
saveRDS(indel.muts.gr, file.path(final.results.dir, "postrank_INDEL_effectorflagged.Rds"))

message("Merging is complete. Final result was written to ", final.results.dir)


# Remove any intermediate files
unlink(intermediate.dir, recursive = TRUE, force = TRUE)
unlink(file.path(opt$`output-dir`, "ranking"), recursive = TRUE, force = TRUE)
file.remove(c(file.path(final.results.dir, "postrank_SNV_effectorflag_byrule.Rds"),
    file.path(final.results.dir, "postrank_INDEL_effectorflag_byrule.Rds")))
