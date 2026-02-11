#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(HDF5Array)
library(ggplot2)

parser <- ArgumentParser(description = "FRASER 2.0 Statistical Analysis")
parser$add_argument("--fds_dir", required = TRUE, help = "Base directory containing the output folder")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--pval_cutoff", type = "double", default = 0.05)
parser$add_argument("--z_cutoff", type = "double")
parser$add_argument("--delta_psi_cutoff", type = "double", default = 0.3)
parser$add_argument("--min_count", type = "integer", default = 5)
parser$add_argument("--nthreads", type = "integer", default = 1)
args <- parser$parse_args()

# --- 1. Parallelization Setup ---
# Use the user-provided nthreads for all BiocParallel operations
bp <- MulticoreParam(workers = args$nthreads)
register(bp)

# --- 2. Load the Dataset ---
options(delayedArray.block.size = 1e9)
fds_name <- paste0("FRASER_", args$cohort_id)
message(paste0("Loading Fraser Data Set: ", fds_name))
fds <- loadFraserDataSet(dir = file.path(args$fds_dir, fds_name), name = fds_name)

# REPAIR: If spliceSiteID is missing, FRASER 2.0 cannot calculate Jaccard/PSI
if(! "spliceSiteID" %in% colnames(mcols(fds, type="ss"))){
  message("Manually injecting Splice Site IDs...")

  # This internal FRASER call generates the mapping without needing the full constructor
  # It populates the 'spliceSiteCoords' slot which calculatePSIValues needs
  fds <- FRASER:::annotateSpliceSite(fds)

  # We also need to ensure the Jaccard-specific metadata is initialized
  # Often, this is what's missing when the merge.data.table fails
  if(is.null(mcols(fds, type="j")$startID)){
    # This maps junctions to the spliceSiteIDs we just generated
    fds <- FRASER:::updateIndices(fds)
  }
}

# Create QC directory
dir.create("qc_plots", showWarnings = FALSE)

# --- 3. Calculate Filter Values FIRST (required for plotting) ---
message("Calculating filter expression values...")
fds <- filterExpressionAndVariability(fds, minDeltaPsi = 0.0, filter = FALSE)

# --- 4. Generate QC Plots ---
# DOWNSAMPLING FOR PLOTS: Use 30,000 random junctions for QC to keep it fast
set.seed(42)
plot_idx <- sample(nrow(fds), min(nrow(fds), 30000))
fds_plot_subset <- fds[plot_idx, ]

message("Generating QC plots using downsampled subset...")
png("qc_plots/filter_expression.png", width = 1200, height = 1200, res = 150)
plotFilterExpression(fds_plot_subset, bins = 100)
dev.off()

# --- 5. Apply Filtering ---
message("Applying filtering based on calculated values...")
fds_filtered <- fds[mcols(fds, type = "j")[, "passed"], ]

# Dimensionality Message
raw_dim <- nrow(fds)
filtered_dim <- nrow(fds_filtered)

message(paste0("\n--- Filtering Summary ---"))
message(paste0("Original junctions:   ", raw_dim))
message(paste0("Filtered junctions:   ", filtered_dim))
message(paste0("Reduction:            ", round((1 - (filtered_dim / raw_dim)) * 100, 2), "%"))

# --- 6. Hyperparameter Optimization ---
# Optimization must run on the filtered set
message("Optimizing encoding dimension (q)...")
opt_q <- bestQ(fds_filtered, type = "jaccard")

png("qc_plots/best_q_optimization.png", width = 1200, height = 1200, res = 150)
plotEncDimSearch(fds_filtered, type = "jaccard")
dev.off()

# --- 7. Fitting ---
message(paste0("Fitting FRASER model with q = ", opt_q, "..."))
fds_fit <- FRASER(fds_filtered, q = opt_q, type = "jaccard", BPPARAM = bp)

# QQ Plot
message("Generating QQ plot...")
png("qc_plots/qq_plot.png", width = 1200, height = 1200, res = 150)
plotQQ(fds_fit, type = "jaccard")
dev.off()

# --- 8. Annotation ---
message("Annotating results with gene information...")
fds_fit <- annotateRangesWithTxDb(fds_fit,
                                  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                  orgDb = org.Hs.eg.db)

# --- 9. Results & Compressed Export ---
message("Extracting results...")
res <- results(fds_fit,
               padjCutoff = args$pval_cutoff,
               deltaPsiCutoff = args$delta_psi_cutoff,
               zScoreCutoff = args$z_cutoff,
               minCount = args$min_count)

# Extract all results for Jaccard
res_all <- results(fds_fit, padjCutoff = 1, deltaPsiCutoff = 0, minCount = 0)

message("Saving results...")
write_csv(as.data.frame(res), paste0(args$cohort_id, ".significant.csv"))
write_csv(as.data.frame(res_all), paste0(args$cohort_id, ".all_results.csv.gz"))

# --- 10. Final Save ---
message("Saving final FRASER object...")
saveFraserDataSet(fds_fit, dir = getwd(), name = paste0(args$cohort_id, "_final"))
message("Analysis Complete.")