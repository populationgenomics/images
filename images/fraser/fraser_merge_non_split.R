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
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)

# 4. Merge Split Counts
# This step updates the fds with the junctions identified in Step 3
message("Merging all split-read counts...")
fds <- getSplitReadCountsForAllSamples(fds, recount = FALSE)

# 5. Merge Non-Split Counts
# We use the filtered ranges from Step 3 to define the 'at-site' junctions
message("Merging all non-split-read counts from cache...")
split_count_ranges <- readRDS(args$filtered_ranges_path)
fds <- getNonSplitReadCountsForAllSamples(
  fds = fds,
  splitCountRanges = split_count_ranges,
  minAnchor = 5,
  recount = FALSE # Crucial: FALSE ensures it uses the .h5 files from the cache
)

# 6. Statistical Normalization (PSI and Jaccard)
# This is the 'FRASER 2.0' logic. It calculates the metrics used for outlier detection.
message("Calculating PSI and Jaccard Index values...")
fds <- calculatePSIValues(fds, BPPARAM = bp)

# 7. Final Save
# This creates the complete 'fitted' dataset that the analysis script will use
message("Saving final merged FraserDataSet...")
fds <- saveFraserDataSet(fds)

message("Merge non-split and PSI calculation complete.")