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
file.copy(args$fds_path, file.path(save_dir, "fds-object.RDS"), overwrite = TRUE)

# 2. Force HDF5 and Configure Parallelism
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = 0) # Changed from -1 to 0 to be safer

# 3. Load and Prune IMMEDIATELY
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_dir_name)

# 4. Update Metadata and Strand
# Important: Ensure the strand matches your actual lab protocol
strandSpecific(fds) <- 0

# SUBSET FIRST: This is the most critical memory-saving step.
# By subsetting here, we drop the metadata of all other samples.
fds <- fds[, fds$sampleID == args$sample_id]
colData(fds)$bamFile <- args$bam_path

# Force Garbage Collection to reclaim memory from the full cohort load
gc()

# 5. Run counting with minimal overhead
# We use recount=FALSE if possible, but keeping TRUE as per your requirement.
# junctionId is set to NULL to ensure a fresh scan of this specific BAM.
fds <- countRNAData(
  fds,
  sampleIds = args$sample_id,
  recount = TRUE,
)

# 6. Verification
# FRASER writes split counts to: {work_dir}/cache/splitCounts/splitCounts-{sample_id}.RDS
expected_out <- file.path(args$work_dir, "cache", "splitCounts", paste0("splitCounts-", args$sample_id, ".RDS"))

if (file.exists(expected_out)) {
    message("Success: ", expected_out)
} else {
    stop("Counting failed to produce output.")
}
