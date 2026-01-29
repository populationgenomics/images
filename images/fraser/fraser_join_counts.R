#!/usr/bin/env Rscript
library(argparse)
library(FRASER)
library(SummarizedExperiment)
library(HDF5Array)

parser <- ArgumentParser(description = "Join Split and Non-Split Counts")
parser$add_argument("--fds_path", required = TRUE)
parser$add_argument("--cohort_id", required = TRUE)
parser$add_argument("--work_dir", default = "/io/work")
parser$add_argument("--filtered_ranges_path", required = TRUE) # Added this back
parser$add_argument("--nthreads", type = "integer", default = 1)
args <- parser$parse_args()

fds_name <- paste0("FRASER_", args$cohort_id)
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
# We need to ensure the se.rds anchor exists for the non-split files moved from Step 5
non_split_count_ranges <- readRDS(args$filtered_ranges_path)

if(!file.exists(file.path(non_split_dir, "se.rds"))){
    message("Creating HDF5 anchor for non-split counts (rawCountsSS)...")
    # This shim tells FRASER how to read the loose .h5 files as a single assay
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
message("Calculating PSI values using both assays...")
# Link the non-split counts now that the anchor is ready
non_split_se <- loadHDF5SummarizedExperiment(dir = non_split_dir)
assay(fds, "rawCountsSS", withDimnames=FALSE) <- assay(non_split_se[, samples(fds)])

# Calculate ratios (PSI)
bp <- MulticoreParam(workers = args$nthreads)
fds <- calculatePSIValues(fds, BPPARAM = bp)

# 5. Save the final integrated object
fds <- saveFraserDataSet(fds)
