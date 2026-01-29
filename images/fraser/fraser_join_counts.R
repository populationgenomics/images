#!/usr/bin/env Rscript
library(argparse)
library(FRASER)
library(SummarizedExperiment)
library(HDF5Array)

parser <- ArgumentParser(description = "Join Split and Non-Split Counts")
parser$add_argument("--fds_path", required = TRUE)
parser$add_argument("--cohort_id", required = TRUE)
parser$add_argument("--filtered_ranges_path", required = TRUE)
parser$add_argument("--work_dir", default = "/io/work")
args <- parser$parse_args()

fds_name <- paste0("FRASER_", args$cohort_id)
# Paths to the HDF5 anchor directories
split_dir <- file.path(args$work_dir, "savedObjects", fds_name, "splitCounts")
non_split_dir <- file.path(args$work_dir, "savedObjects", fds_name, "nonSplitCounts")

# 1. Load the FDS object
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)

# 2. Re-anchor Split Counts (The 'j' counts)
# If the splitCounts folder was moved/localized, fds might lose the pointer.
# We reload the HDF5 SummarizedExperiment and inject it back into fds.
if(file.exists(file.path(split_dir, "se.rds"))){
    message("Re-anchoring split counts (rawCountsJ)...")
    split_se <- loadHDF5SummarizedExperiment(dir = split_dir)
    assay(fds, "rawCountsJ", withDimnames=FALSE) <- assay(split_se)
} else {
    stop("CRITICAL ERROR: splitCounts/se.rds not found. Cannot calculate PSI.")
}

# 3. Handle Non-Split Counts (The 'ss' counts)
# Use the shim logic from before to ensure the 'rawCountsSS' is present
non_split_count_ranges <- readRDS(args$filtered_ranges_path)
if(!file.exists(file.path(non_split_dir, "se.rds"))){
    message("Creating anchor for non-split counts (rawCountsSS)...")
    anchor_se <- SummarizedExperiment(
        assays = list(rawCountsSS = matrix(0, nrow=length(non_split_count_ranges),
                                           ncol=length(samples(fds)),
                                           dimnames=list(NULL, samples(fds)))),
        rowRanges = non_split_count_ranges
    )
    saveHDF5SummarizedExperiment(anchor_se, dir=non_split_dir, replace=TRUE, prefix="")
}

# 4. Calculate PSI
# Now that both 'rawCountsJ' and 'rawCountsSS' are linked, this will succeed
message("Calculating PSI values...")
fds <- calculatePSIValues(fds)

message("Saving integrated FraserDataSet...")
fds <- saveFraserDataSet(fds)
