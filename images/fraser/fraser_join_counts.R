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

# Paths to the HDF5 anchor directories
split_dir <- file.path(saveDir, "splitCounts")
non_split_dir <- file.path(saveDir, "nonSplitCounts")

# Paths to the RDS files generated in previous steps
gRangesSplitPath <- file.path(args$work_dir, "g_ranges_split_counts.RDS")
spliceSitePath <- file.path(args$work_dir, "g_ranges_non_split_counts.RDS")

# 1. Load the FDS object
message("Loading Fraser Data Set...")
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)

# Load the ranges
splitCounts_gRanges <- readRDS(gRangesSplitPath)
spliceSiteCoords <- readRDS(spliceSitePath)

# 2. Get splitReads counts
# The HDF5 file is named 'rawCountsJ.h5' inside the splitCounts subdir
splitCounts_h5_path <- file.path(split_dir, "rawCountsJ.h5")
message(paste("Loading split counts from:", splitCounts_h5_path))

splitCounts_h5 <- HDF5Array::HDF5Array(splitCounts_h5_path, "rawCountsJ")
splitCounts_se <- SummarizedExperiment(
  colData = colData(fds),
  rowRanges = splitCounts_gRanges,
  assays = list(rawCountsJ = splitCounts_h5)
)

# 3. Get nonSplitRead counts
# The HDF5 file is named 'rawCountsSS.h5' inside the nonSplitCounts subdir
nonSplitCounts_h5_path <- file.path(non_split_dir, "rawCountsSS.h5")
message(paste("Loading non-split counts from:", nonSplitCounts_h5_path))

nonSplitCounts_h5 <- HDF5Array::HDF5Array(nonSplitCounts_h5_path, "rawCountsSS")
nonSplitCounts_se <- SummarizedExperiment(
  colData = colData(fds),
  rowRanges = spliceSiteCoords,
  assays = list(rawCountsSS = nonSplitCounts_h5)
)

# 4. Add Counts to FRASER object
# This function automatically populates the internal slots correctly
message("Joining assays into FDS object...")
fds <- addCountsToFraserDataSet(
  fds = fds,
  splitCounts = splitCounts_se,
  nonSplitCounts = nonSplitCounts_se
)

# 5. Calculate PSI values
message("Calculating PSI values...")
bp <- MulticoreParam(workers = args$nthreads)
fds <- calculatePSIValues(fds, BPPARAM = bp)

# 6. Save final FRASER object
message("Saving integrated FDS...")
fds <- saveFraserDataSet(fds)

message("FRASER join complete.")
