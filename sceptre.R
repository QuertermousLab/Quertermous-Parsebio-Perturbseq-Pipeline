#!/usr/bin/env Rscript

# Load dependencies
suppressMessages(library(sceptre))
suppressMessages(library(cowplot))
suppressMessages(library(optparse))

# Define command line options
option_list <- list(
  make_option(c("--gene_mat"), type="character", help="Path to gene count matrix (.mtx)"),
  make_option(c("--grna_mat"), type="character", help="Path to gRNA count matrix (.mtx)"),
  make_option(c("--all_genes"), type="character", help="Path to all_genes.csv"),
  make_option(c("--all_grnas"), type="character", help="Path to all_guides.csv"),
  make_option(c("--grna_target_df"), type="character", help="Path to grna_target_data_frame.csv"),
  make_option(c("--pos_pairs"), type="character", help="Path to positive control pairs CSV"),
  make_option(c("--output_dir"), type="character", help="Directory to save outputs"),
  make_option(c("--rds_file"), type="character", help="Path to save sceptre RDS object"),
  make_option(c("--n_processors"), type="integer", default=1, help="Number of processors for parallel operations")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Create output directories if needed
if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
if (!dir.exists(dirname(opt$rds_file))) dir.create(dirname(opt$rds_file), recursive = TRUE)

# Import data
grna_target_data_frame <- read.csv(opt$grna_target_df)
sceptre_object <- import_data_from_parse(
  gene_mat_fp = opt$gene_mat,
  grna_mat_fp = opt$grna_mat,
  all_genes_fp = opt$all_genes,
  all_grnas_fp = opt$all_grnas,
  moi = "low",
  grna_target_data_frame = grna_target_data_frame
)

# Construct analysis pairs
positive_control_pairs <- read.csv(opt$pos_pairs)
discovery_pairs_lowmoi <- construct_trans_pairs(
  sceptre_object = sceptre_object,
  positive_control_pairs = positive_control_pairs,
  pairs_to_exclude = "none"
)

# Set analysis parameters
sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs_lowmoi,
  positive_control_pairs = positive_control_pairs,
  side = 'both'
)

# Assign gRNAs
sceptre_object_mixture <- assign_grnas(
  sceptre_object = sceptre_object, 
  method = "mixture",
  parallel = TRUE, 
  n_processors = opt$n_processors
)

# Run QC
sceptre_object <- run_qc(sceptre_object, p_mito_threshold = 0.075)

# Calibration check
sceptre_object <- run_calibration_check(
  sceptre_object, 
  parallel = TRUE, 
  n_processors = opt$n_processors
)
calibration_result <- get_result(sceptre_object, analysis = "run_calibration_check")
write.table(calibration_result, file=file.path(opt$output_dir, "calibration_result.txt"), sep='\t', quote=FALSE)

# Power check
sceptre_object <- run_power_check(
  sceptre_object, 
  parallel = TRUE, 
  n_processors = opt$n_processors
)

# Discovery analysis
sceptre_object <- run_discovery_analysis(
  sceptre_object, 
  parallel = TRUE, 
  n_processors = opt$n_processors
)

# Write results
write_outputs_to_directory(sceptre_object, directory = opt$output_dir)
result <- get_result(sceptre_object, analysis = "run_discovery_analysis")
result$p.adj <- p.adjust(result$p_value, method = 'BH')
write.table(result, file=file.path(opt$output_dir, "discovery_result.txt"), sep='\t', quote=FALSE, row.names=FALSE)

# Save RDS
saveRDS(sceptre_object, file = opt$rds_file)
