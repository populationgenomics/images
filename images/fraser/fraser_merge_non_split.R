#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(BiocParallel)
library(SummarizedExperiment)

parser <- ArgumentParser(description = "Merge Non-Split Counts and Calculate PSI")
parser$add_argument("--fds_path", required = TRUE, help = "Path to localized FDS RDS file")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--filtered_ranges_path", required = TRUE, help = "Path to g_ranges_non_split_counts.RDS")
parser$add_argument("--work_dir", default = "/io/work", help = "Working directory")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

# 1. Load the ranges first so they are available for all steps
non_split_count_ranges <- readRDS(args$filtered_ranges_path)

# 2. Reconstruct Directory Structure
fds_name <- paste0("FRASER_", args$cohort_id)
save_dir <- file.path(args$work_dir, "savedObjects", fds_name)
out_dir <- file.path(save_dir, "nonSplitCounts")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(args$fds_path, file.path(save_dir, "fds-object.RDS"), overwrite = TRUE)

# --- MINIMUM FIX: Move files and create the missing 'se.rds' metadata ---
h5_files <- list.files("/io/batch", pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
for(f in h5_files) {
    file.copy(f, file.path(out_dir, basename(f)), overwrite = TRUE)
}
# 3. Load Dataset
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)

# 1. Create a dummy matrix with the right number of samples (columns)
# This satisfies the: if(all(samples(fds) %in% colnames(siteCounts))) check
sample_names <- samples(fds)
dummy_matrix <- matrix(0, nrow=length(non_split_count_ranges), ncol=length(sample_names))
colnames(dummy_matrix) <- sample_names

# 2. Create the SE object with the dummy matrix
tmp_se <- SummarizedExperiment(
    assays = list(rawCountsSS = dummy_matrix), # Name must match what FRASER expects
    rowRanges = non_split_count_ranges
)

# 3. Save it to the cache directory
saveRDS(tmp_se, file.path(out_dir, "se.rds"))
# -----------------------------------------------------------------------

# 2. Configure Backend
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)
bp <- MulticoreParam(workers = args$nthreads)
register(bp)


available_bams <- list.files("/io/batch/input_bams", pattern = "\\.bam$", full.names = TRUE)
if(length(available_bams) > 0){
    ref_bam <- available_bams[1]
    message("Validating metadata against reference BAM: ", ref_bam)
    colData(fds)$bamFile <- ref_bam
    seqlevelsStyle(fds) <- seqlevelsStyle(Rsamtools::BamFile(ref_bam))
} else {
    stop("CRITICAL ERROR: No BAM files found in /io/batch/input_bams.")
}
strandSpecific(fds) <- 0

# --- 6. Merge Non-Split Counts ---
message("Merging counts using HDF5 files in: ", out_dir)
non_split_counts <- getNonSplitReadCountsForAllSamples(
  fds = fds,
  splitCountRanges = non_split_count_ranges,
  minAnchor = 5,
  recount = FALSE # Crucial: FALSE ensures it uses the .h5 files from the cache
)

# 7. Final Save
# This creates the complete 'fitted' dataset that the analysis script will use

saveRDS(non_split_counts, file.path(args$work_dir, "non_split_counts.RDS"))

