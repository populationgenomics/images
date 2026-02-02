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

# 1. Load the merged FDS
# Python extracted the tar into args$fds_dir.
# Inside is 'savedObjects/FRASER_{cohort_id}/...'
options(delayedArray.block.size = 1e9) # 1GB blocks
fds_name <- paste0("FRASER_", args$cohort_id)
message(paste0("Loading Fraser Data Set: ", fds_name, " from ", args$fds_dir))
fds <- loadFraserDataSet(dir = file.path(args$fds_dir, fds_name), name = fds_name)


# --- 3. Filtering ---
# It is critical to filter before fitting to reduce the size of the latent space matrices
fds <- calculatePSIValues(fds, types = "jaccard", BPPARAM = bp)

# Create QC directory
dir.create("qc_plots", showWarnings = FALSE)

# DOWNSAMPLING FOR PLOTS: Use 30,000 random junctions for QC to keep it fast
set.seed(42)
plot_idx <- sample(nrow(fds), min(nrow(fds), 30000))
fds_plot_subset <- fds[plot_idx, ]

message("Generating QC plots using downsampled subset...")
png("qc_plots/filter_expression.png", width = 1200, height = 1200, res = 150)
plotFilterExpression(fds_plot_subset, bins = 100)
dev.off()

# Apply actual filtering to full dataset
fds <- filterExpressionAndVariability(fds, minDeltaPsi = 0.0, filter = TRUE)

# --- 4. Dimensionality Message ---
raw_dim <- nrow(fds)
fds_filtered <- fds[mcols(fds, type = "j")[, "passed"], ]
filtered_dim <- nrow(fds_filtered)

message(paste0("\n--- Filtering Summary ---"))
message(paste0("Original junctions:   ", raw_dim))
message(paste0("Filtered junctions:   ", filtered_dim))
message(paste0("Reduction:            ", round((1 - (filtered_dim / raw_dim)) * 100, 2), "%"))

# --- 5. Hyperparameter Optimization ---
# Optimization must run on the filtered set
opt_q <- bestQ(fds_filtered, type = "jaccard")

png("qc_plots/best_q_optimization.png", width = 1200, height = 1200, res = 150)
plotEncDimSearch(fds_filtered, type = "jaccard")
dev.off()

# --- 6. Fitting ---
fds_fit <- FRASER(fds_filtered, q = opt_q, type = "jaccard", BPPARAM = bp)

# QQ Plot (also uses a subset internally in FRASER, but we'll be explicit)
png("qc_plots/qq_plot.png", width = 1200, height = 1200, res = 150)
plotQQ(fds_fit, type = "jaccard")
dev.off()

# --- 7. Annotation ---
fds_fit <- annotateRangesWithTxDb(fds_fit,
                                  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                  orgDb = org.Hs.eg.db)

# --- 8. Results & Compressed Export ---
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

# --- 9. Final Save ---
saveFraserDataSet(fds_fit, dir = getwd(), name = paste0(args$cohort_id, "_final"))
message("Analysis Complete.")