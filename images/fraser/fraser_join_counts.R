#!/usr/bin/env Rscript
library(argparse)
library(FRASER)
library(SummarizedExperiment)
library(HDF5Array)
library(BiocParallel)

parser <- ArgumentParser(description = "Join Split and Non-Split Counts")
parser$add_argument("--fds_path", required = TRUE)
parser$add_argument("--cohort_id", required = TRUE)
parser$add_argument("--work_dir", default = "/io/work")
parser$add_argument("--nthreads", type = "integer", default = 1)
args <- parser$parse_args()

fds_name <- paste0("FRASER_", args$cohort_id)
# The orchestrator puts the non-split ranges here:
filtered_ranges_path <- file.path(args$work_dir, "g_ranges_non_split_counts.RDS")

# Paths to the HDF5 anchor directories
split_dir <- file.path(args$work_dir, "savedObjects", fds_name, "splitCounts")
non_split_dir <- file.path(args$work_dir, "savedObjects", fds_name, "nonSplitCounts")

# 1. Load the FDS object
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)

# 2. Re-attach Split Counts (Junctions)
# This prevents the "Missing rawCounts for 'j'" error
if(file.exists(file.path(split_dir, "se.rds"))){
    message("Linking split counts (rawCountsJ)...")
    split_se <- loadHDF5SummarizedExperiment(dir = split_dir)
    assay(fds, "rawCountsJ", withDimnames=FALSE) <- assay(split_se[, samples(fds)])
}

# 3. Handle Non-Split Counts (Splice Sites - 'ss')
# We must use the ranges generated in the 'Merge Split' step
if(!file.exists(filtered_ranges_path)){
    stop(paste("Required non-split ranges not found at:", filtered_ranges_path))
}
non_split_count_ranges <- readRDS(filtered_ranges_path)

if(!file.exists(file.path(non_split_dir, "se.rds"))){
    message("Creating HDF5 anchor for non-split counts (rawCountsSS)...")

    # Create the container matching the dimensions FRASER expects
    anchor_se <- SummarizedExperiment(
        assays = list(rawCountsSS = matrix(0,
                                           nrow = length(non_split_count_ranges),
                                           ncol = length(samples(fds)),
                                           dimnames = list(NULL, samples(fds)))),
        rowRanges = non_split_count_ranges
    )
    # This creates the se.rds that FRASER needs to 'see' the .h5 files
    saveHDF5SummarizedExperiment(anchor_se, dir = non_split_dir, replace = TRUE, prefix = "")
}

# 4. Final Join and PSI Calculation
message("Loading non-split counts...")
non_split_se <- loadHDF5SummarizedExperiment(dir = non_split_dir)
assay(fds, "rawCountsSS", withDimnames=FALSE) <- assay(non_split_se[, samples(fds)])

# Calculate ratios (PSI)
message("Calculating PSI values...")
bp <- MulticoreParam(workers = args$nthreads)
fds <- calculatePSIValues(fds, BPPARAM = bp)

# 5. Save the final integrated object
fds <- saveFraserDataSet(fds)
