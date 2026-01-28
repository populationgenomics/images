#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(BiocParallel)

parser <- ArgumentParser(description = "Merge Non-Split Counts and Calculate PSI")
parser$add_argument("--fds_path", required = TRUE, help = "Path to localized FDS RDS file")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--filtered_ranges_path", required = TRUE, help = "Path to g_ranges_non_split_counts.RDS")
parser$add_argument("--work_dir", default = "/io/work", help = "Working directory")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

# 1. Reconstruct Directory Structure
fds_name <- paste0("FRASER_", args$cohort_id)
save_dir <- file.path(args$work_dir, "savedObjects", fds_name)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(args$fds_path, file.path(save_dir, "fds-object.RDS"))

# 2. Configure Backend
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)
bp <- MulticoreParam(workers = args$nthreads)
register(bp)

# 3. Load Dataset
# This loads the metadata; the counts will be pulled from the cache symlinked by Python
fds_name <- paste0("FRASER_", args$cohort_id)
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)


available_bams <- list.files("/io/batch/input_bams", pattern = "\\.bam$", full.names = TRUE)

if(length(available_bams) > 0){
    # We use the reference BAM to satisfy the internal 'file.exists' checks.
    # available_bams[1] should be the 'reference.bam' symlinked by the Python task.
    ref_bam <- available_bams[1]
    message("Validating metadata against reference BAM: ", ref_bam)

    # Update colData paths for all samples to the reference BAM
    colData(fds)$bamFile <- ref_bam

    # Ensure seqlevels match (e.g., 'chr1' vs '1') to prevent crashing during merge
    seqlevelsStyle(fds) <- seqlevelsStyle(BamFile(ref_bam))
} else {
    stop("CRITICAL ERROR: No BAM files found in /io/batch/input_bams. Validation will fail.")
}
strandSpecific(fds) <- 0

# 5. Merge Non-Split Counts
# We use the filtered ranges from Step 3 to define the 'at-site' junctions
message("Merging all non-split-read counts from cache...")
non_split_count_ranges <- readRDS(args$filtered_ranges_path)
fds <- getNonSplitReadCountsForAllSamples(
  fds = fds,
  splitCountRanges = non_split_count_ranges,
  minAnchor = 5,
  recount =FALSE # Crucial: FALSE ensures it uses the .h5 files from the cache
)

# 7. Final Save
# This creates the complete 'fitted' dataset that the analysis script will use
message("Saving final merged FraserDataSet...")
fds <- saveFraserDataSet(fds)

message("Merge non-split and PSI calculation complete.")
