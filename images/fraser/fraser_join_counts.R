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

# Set number of threads
register(MulticoreParam(args$nthreads))

fds_name <- paste0("FRASER_", args$cohort_id)

# Load the FDS object
message(paste("Loading FraserDataSet from:", args$work_dir))
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)

# calculatePSIValues will look for HDF5 arrays in:
# 1. work_dir/savedObjects/fds_name/splitCounts
# 2. work_dir/savedObjects/fds_name/nonSplitCounts
message("Calculating PSI values using HDF5-backed assays...")
fds <- calculatePSIValues(fds)

message("Saving integrated FraserDataSet...")
fds <- saveFraserDataSet(fds)
