#!/usr/bin/env Rscript

# Set memory limit at the environment level before loading heavy libraries
# 11GB limit
Sys.setenv("R_MAX_VSIZE" = "14Gb")

library(argparse)
library(FRASER)
library(BiocParallel)

parser <- ArgumentParser(description = "Count Split Reads for a Single Sample")
parser$add_argument("--fds_path", required = TRUE, help = "Path to localized FDS RDS file")
parser$add_argument("--bam_path", required = TRUE, help = "Path to the specific BAM file for this sample")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--sample_id", required = TRUE, help = "Sample ID to count")
parser$add_argument("--work_dir", default = "/io/work", help = "Working directory")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

# 1. Reconstruct the FRASER directory structure so loadFraserDataSet works
# FRASER expects: {work_dir}/savedObjects/FRASER_{cohort_id}/fds-object.RDS
fds_dir_name <- paste0("FRASER_", args$cohort_id)
save_dir <- file.path(args$work_dir, "savedObjects", fds_dir_name)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(args$fds_path, file.path(save_dir, "fds-object.RDS"))

# 2. Force HDF5 and Configure Parallelism
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)

# Configure MulticoreParam with the 11GB limit in mind
# Note: nthreads shares the 11GB pool
bpparam <- MulticoreParam(
    workers = args$nthreads,
    stop.on.error = TRUE
)
register(bpparam)

# 3. Load the dataset
# Note: name must match what was used in fraser_init.R
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_dir_name)

# 4. Update Metadata and Strand
# Important: Ensure the strand matches your actual lab protocol
strandSpecific(fds) <- 0

# Subset the FDS object to ONLY this sample.
# This prevents the worker from erroring out when looking for other samples' BAMs.
if (!(args$sample_id %in% colData(fds)$sampleID)) {
    stop(paste("Sample", args$sample_id, "not found in the provided FDS object."))
}
fds <- fds[, fds$sampleID == args$sample_id]

# Update the localized BAM path in colData
colData(fds)$bamFile <- args$bam_path

# 5. Run counting
# We pass BPPARAM explicitly to ensure the worker uses the allocated threads
fds <- countRNAData(
  fds,
  sampleIds = args$sample_id,
  recount = TRUE,
  BPPARAM = bpparam
)

# 6. Verification
# FRASER writes split counts to: {work_dir}/cache/splitCounts/splitCounts-{sample_id}.RDS
expected_out <- file.path(args$work_dir, "cache", "splitCounts", paste0("splitCounts-", args$sample_id, ".RDS"))

if (file.exists(expected_out)) {
    message("Successfully created split counts at: ", expected_out)
} else {
    # List files in the directory to help debug if it fails
    found_files <- list.files(file.path(args$work_dir, "cache", "splitCounts"))
    stop(paste0(
        "Split counts file not found. Expected: ", expected_out,
        "\nFound in cache: ", paste(found_files, collapse=", ")
    ))
}
