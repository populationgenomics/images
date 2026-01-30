#!/usr/bin/env Rscript
library(argparse)
library(FRASER)
library(SummarizedExperiment)
library(HDF5Array)

parser <- ArgumentParser(description = "Join Split and Non-Split Counts")
parser$add_argument("--fds_path", required = TRUE)
parser$add_argument("--cohort_id", required = TRUE)
parser$add_argument("--work_dir", default = "/io/work")
parser$add_argument("--nthreads", type = "integer", default = 1)
args <- parser$parse_args()

fds_name <- paste0("FRASER_", args$cohort_id)
saveDir <- file.path(args$work_dir, "savedObjects", fds_name)

# 1. Load the FDS object
message("Loading Fraser Data Set from RDS...")
fds <- readRDS(args$fds_path)

# Update the internal directory path to match the current work directory
workingDir(fds) <- saveDir

# 2. Load split counts
message("Loading split counts...")
split_counts_path <- file.path(args$work_dir, "g_ranges_split_counts.RDS")
if(!file.exists(split_counts_path)){
    stop(paste("Missing split counts at:", split_counts_path))
}
splitCounts_se <- readRDS(split_counts_path)

# CRITICAL: Materialize the assays if they're HDF5-backed with broken paths
message("Materializing split counts assays...")
for(assay_name in names(assays(splitCounts_se))){
    assays(splitCounts_se)[[assay_name]] <- as.matrix(assays(splitCounts_se)[[assay_name]])
}

# 3. Load merged non-split counts
message("Loading merged non-split counts...")
merged_non_split_dir <- file.path(args$work_dir, "merged_non_split_counts")
if(!dir.exists(merged_non_split_dir)){
    stop(paste("Missing merged non-split counts directory at:", merged_non_split_dir))
}
nonSplitCounts_se <- loadHDF5SummarizedExperiment(dir = merged_non_split_dir)

# CRITICAL: Materialize non-split counts assays too
message("Materializing non-split counts assays...")
for(assay_name in names(assays(nonSplitCounts_se))){
    assays(nonSplitCounts_se)[[assay_name]] <- as.matrix(assays(nonSplitCounts_se)[[assay_name]])
}

# 4. Add counts to FRASER object
message("Joining assays into FDS object...")
fds <- addCountsToFraserDataSet(
  fds = fds,
  splitCounts = splitCounts_se,
  nonSplitCounts = nonSplitCounts_se
)

# 5. Save final FRASER object
message("Saving final integrated FDS...")
fds <- saveFraserDataSet(fds)

message("FRASER join complete.")
