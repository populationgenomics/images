#!/usr/bin/env Rscript

# Set memory limit - increased slightly to allow for HDF5 overhead
Sys.setenv("R_MAX_VSIZE" = "16Gb")

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

# Create the cache directory WITHOUT the fds_name subdirectory
# This matches what the merge script expects
cache_dir <- file.path(args$work_dir, "cache", "nonSplicedCounts")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Configure Parallelism and HDF5
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)
# Use MulticoreParam if threads > 1, else Serial
if(args$nthreads > 1){
    bpparam <- MulticoreParam(workers = args$nthreads)
} else {
    bpparam <- SerialParam()
}

# 3. Load Dataset and Coordinates
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)
filtered_coords <- readRDS(args$coords_path) # This is your splice_site_coords.RDS
# Add this before calling countNonSplicedReads
colData(fds)$bamFile <- args$bam_path

# Set strand specificity
strandSpecific(fds) <- 2

# 5. Run Non-Split Counting
# This writes the .h5 or .RDS file into the cache_dir created above
message(paste("Counting non-split reads for sample:", args$sample_id))

sample_result <- countNonSplicedReads(args$sample_id,
                                      splitCountRanges = NULL,
                                      fds = fds,
                                      NcpuPerSample = args$nthreads,
                                      minAnchor=5,
                                      recount= TRUE,
                                      spliceSiteCoords=filtered_coords)
# 6. Verification
# The H5 file should now be at cache_dir/nonSplicedCounts-{sample_id}.h5
expected_h5  <- file.path(cache_dir, paste0("nonSplicedCounts-", args$sample_id, ".h5"))
expected_rds <- file.path(cache_dir, paste0("nonSplicedCounts-", args$sample_id, ".RDS"))

if (file.exists(expected_h5)) {
    message("Success: Created non-split counts (HDF5) at ", expected_h5)
} else if (file.exists(expected_rds)) {
    message("Success: Created non-split counts (RDS) at ", expected_rds)
} else {
    # Debugging: List files in the directory to see what actually happened
    if(dir.exists(cache_dir)){
        message("Files found in cache dir: ", paste(list.files(cache_dir), collapse=", "))
    }

    if(!file.exists(paste0(args$bam_path, ".bai"))){
        stop("BAM Index (.bai) missing. FRASER cannot perform random access.")
    }
    stop(paste("Counting failed. Output not found for sample:", args$sample_id))
}