#!/usr/bin/env Rscript

# Set memory limit - increased slightly to allow for HDF5 overhead
Sys.setenv("R_MAX_VSIZE" = "16Gb")

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
options("FRASER.maxJunctionsNoHDF5" = -1)
# Use MulticoreParam if threads > 1, else Serial
if(args$nthreads > 1){
    bpparam <- MulticoreParam(workers = args$nthreads)
} else {
    bpparam <- SerialParam()
}

# 3. Load and Prune IMMEDIATELY
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_dir_name)

strandSpecific(fds) <- 2

# SUBSET FIRST: This is the most critical memory-saving step.
# By subsetting here, we drop the metadata of all other samples.
fds <- fds[, fds$sampleID == args$sample_id]

# Validate the BAM path - Ensure the R script sees what Hail localized
if(!file.exists(args$bam_path)){
    stop(paste("BAM file not found at:", args$bam_path))
}
colData(fds)$bamFile <- args$bam_path

# 4. Count Split Reads
message(paste("Starting split read counting for sample:", args$sample_id))

# In FRASER 2.0, we use getSplitReadCountsForAllSamples with recount=TRUE.
# This writes an RDS file to the cache which we will harvest in the merge step.
fds <- getSplitReadCountsForAllSamples(
    fds,
    recount = TRUE,
    BPPARAM = bpparam
)

# 5. Verification
# FRASER saves individual counts to: cache/splitCounts/splitCounts-{sample_id}.RDS
expected_out <- file.path(args$work_dir, "cache", "splitCounts", paste0("splitCounts-", args$sample_id, ".RDS"))

if (file.exists(expected_out)) {
    message("Success: Created split counts at ", expected_out)
} else {
    # Check for common Bioinformatic failures
    if(!file.exists(paste0(args$bam_path, ".bai"))){
        stop("BAM Index (.bai) missing. FRASER cannot perform random access counting.")
    }
    stop("Counting failed. The RDS file was not generated in the cache.")
}
