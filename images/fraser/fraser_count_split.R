#!/usr/bin/env Rscript

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
register(MulticoreParam(workers = args$nthreads))

# 3. Load the dataset
# Note: name must match what was used in fraser_init.R
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_dir_name)

# 4. Update Metadata and Strand
# Important: Ensure the strand matches your actual lab protocol
strandSpecific(fds) <- 0

# Use colData to update the localized BAM path
if (args$sample_id %in% colData(fds)$sampleID) {
    colData(fds)[colData(fds)$sampleID == args$sample_id, "bamFile"] <- args$bam_path
} else {
    stop("Sample ID not found in fds object.")
}
# 5. Run counting for the specific sample
# countRNAData with sampleId filter is the standard way to run parallel counts
fds <- countRNAData(
  fds,
  sampleIds = args$sample_id,
  NcpuPerSample = args$nthreads,
  recount = TRUE
)

# 6. Verification
# FRASER writes split counts to: {work_dir}/cache/splitCounts/splitCounts-{sample_id}.RDS
expected_out <- file.path(args$work_dir, "cache", "splitCounts", paste0("splitCounts-", args$sample_id, ".RDS"))

if (file.exists(expected_out)) {
    message("Successfully created split counts at: ", expected_out)
} else {
    stop("Split counts file was not found at expected location: ", expected_out)
}
