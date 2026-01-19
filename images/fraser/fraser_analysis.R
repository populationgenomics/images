#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
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

# Force HDF5 for large matrix operations
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)

z_cutoff <- if (is.null(args$z_cutoff)) NA else args$z_cutoff
bp <- MulticoreParam(workers = args$nthreads)
register(bp)

# 1. Load the merged FDS
# Python extracted the tar into args$fds_dir.
# Inside is 'savedObjects/FRASER_{cohort_id}/...'
fds_name <- paste0("FRASER_", args$cohort_id)
fds <- loadFraserDataSet(dir = args$fds_dir, name = fds_name)

# 2. Filter Expression and Variability
# This reduces the number of tests and improves power
message("Filtering junctions...")
fds <- filterExpressionAndVariability(fds, minDeltaPsi = 0.0, filter = TRUE)

dir.create("plots/misc", recursive = TRUE, showWarnings = FALSE)
png("plots/misc/filter_expression.png", width = 2000, height = 2000, res = 300)
plotFilterExpression(fds, bins = 100)
dev.off()

# 3. Hyperparameter Optimization (bestQ) and Fitting
# FRASER uses an autoencoder to control for latent confounders (batch effects)
message("Finding optimal hyperparameter Q and fitting autoencoder...")
# In FRASER 2.0, we typically fit psi5, psi3, and jaccard
fds <- FRASER(fds, q = c(psi5=NA, psi3=NA, jaccard=NA), BPPARAM = bp)

# 4. Annotation
message("Annotating results with UCSC hg38...")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgDb <- org.Hs.eg.db
fds <- annotateRangesWithTxDb(fds, txdb = txdb, orgDb = orgDb)

# 5. Extract Results
message("Exporting results...")
res <- results(fds,
               padjCutoff = args$pval_cutoff,
               deltaPsiCutoff = args$delta_psi_cutoff,
               zScoreCutoff = z_cutoff,
               minCount = args$min_count)

res_all <- results(fds, padjCutoff = 1, deltaPsiCutoff = 0, minCount = 0)

write_csv(as.data.frame(res), "results.significant.csv")
write_csv(as.data.frame(res_all), "results.all.csv")

# 6. Generate Diagnostic Plots
message("Generating Volcano and Q-Q plots...")
dir.create("plots/volcano", recursive = TRUE)
for(type in c("psi5", "psi3", "jaccard")){
    png(paste0("plots/volcano/volcano_", type, ".png"), width = 2000, height = 2000, res = 300)
    print(plotVolcano(fds, type = type, sampleID = sampleIDs(fds)[1])) # Plots first sample as example
    dev.off()
}

# 7. Final Save
# We save this in a generic folder so Python can easily grab it
saveFraserDataSet(fds, dir = "output", name = paste0("Final_", args$cohort_id))

# Summary statistics for the log
sink("statistics_summary.txt")
cat("Cohort ID:", args$cohort_id, "\n")
cat("Total Samples:", length(sampleIDs(fds)), "\n")
cat("Significant Events (adj P <", args$pval_cutoff, "):", nrow(res), "\n")
sink()

message("Analysis complete.")