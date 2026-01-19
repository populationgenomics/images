#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(BiocParallel)

parser <- ArgumentParser(description = "Count Non-Split Reads for a Single Sample")
parser$add_argument("--fds_path", required = TRUE, help = "Path to localized FDS RDS file")
parser$add_argument("--bam_path", required = TRUE, help = "Path to the localized BAM file")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--sample_id", required = TRUE, help = "Sample ID")
parser$add_argument("--coords_path", required = TRUE, help = "Path to splice_site_coords.RDS")
parser$add_argument("--work_dir", default = "/io/work", help = "Working directory")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

# 1. Reconstruct Directory Structure
fds_name <- paste0("FRASER_", args$cohort_id)
save_dir <- file.path(args$work_dir, "savedObjects", fds_name)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(args$fds_path, file.path(save_dir, "fds-object.RDS"))

# Create the specific cache directory for non-split counts
cache_dir <- file.path(args$work_dir, "cache", "nonSplicedCounts", fds_name)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Configure Parallelism and HDF5
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)
register(MulticoreParam(workers = args$nthreads))

# 3. Load Dataset and Coordinates
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)
splice_site_coords <- readRDS(args$coords_path)

# 4. Synchronize BAM Path
# Ensure the FDS object points to the current local BAM location on the worker
bamData(fds)[bamData(fds)$sampleID == args$sample_id, "bamFile"] <- args$bam_path

# 5. Run Non-Split Counting
# This writes the .h5 or .RDS file into the cache_dir created above
message(paste("Counting non-split reads for sample:", args$sample_id))
countNonSplicedReads(
  sampleID = args$sample_id,
  fds = fds,
  spliceSiteCoords = splice_site_coords,
  splitCountRanges = NULL, # Already handled in merge_split
  minAnchor = 5,
  NcpuPerSample = args$nthreads,
  recount = TRUE
)

# 6. Verification
# The Python script uses a 'find' command to move this output.
# We just need to check if the file was created in the cache.
out_file <- list.files(cache_dir, pattern = paste0("nonSplicedCounts-", args$sample_id))
if (length(out_file) > 0) {
    message("Successfully created non-split counts: ", out_file[1])
} else {
    stop("Non-split counts file was not found in: ", cache_dir)
}