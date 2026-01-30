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

# Paths to the RDS files generated in previous steps
gRangesSplitPath <- file.path(args$work_dir, "g_ranges_split_counts.RDS")
spliceSitePath <- file.path(args$work_dir, "splice_site_coords.RDS")

# 1. Load the FDS object DIRECTLY from RDS (not using loadFraserDataSet)
message("Loading Fraser Data Set from RDS...")
fds <- readRDS(args$fds_path)

# Update the internal directory path to match the current work directory
workingDir(fds) <- saveDir

# Load the ranges
splitCounts_gRanges <- readRDS(gRangesSplitPath)
spliceSiteCoords <- readRDS(spliceSitePath)

# 2. Get splitReads counts
message("Loading split counts...")
split_se_path <- file.path(saveDir, "splitCounts", "se.rds")
if(!file.exists(split_se_path)){
    stop(paste("Missing splitCounts anchor at:", split_se_path))
}
splitCounts_se <- readRDS(split_se_path)

# 3. Get nonSplitRead counts
message("Loading non-split counts...")
non_split_se_path <- file.path(saveDir, "nonSplitCounts", "se.rds")
if(!file.exists(non_split_se_path)){
    stop(paste("Missing nonSplitCounts anchor at:", non_split_se_path))
}
nonSplitCounts_se <- readRDS(non_split_se_path)

# 4. Add Counts to FRASER object
# This function automatically populates the internal slots correctly
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
