library(argparse)
library(FRASER)
library(BiocParallel)

parser <- ArgumentParser(description = "Initialize FRASER Data Set")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--sample_map", required = TRUE, help = "Path to CSV with sample_id and bam_path")
parser$add_argument("--work_dir", default = "/io/work", help = "Working directory for FRASER")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

# 1. Read the static mapping created by the Python job
sample_map <- read.csv(args$sample_map, stringsAsFactors = FALSE)
# Define the dataset name consistently with the Python expectations
fds_name <- paste0("FRASER_", args$cohort_id)


sample_table <- DataFrame(
  sampleID  = as.character(sample_map$sample_id),
  bamFile   = as.character(sample_map$bam_path),
  pairedEnd = TRUE
)
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)

fds <- FraserDataSet(
  colData    = sample_table,
  workingDir = args$work_dir,
  name       = fds_name
)

#Setup parallel execution
bp <- MulticoreParam(workers = args$nthreads)
register(bp)


# We calculate initial metadata/PSI skeleton

fds <- saveFraserDataSet(fds)

#Print location for Python to capture if needed
fds_save_path <- file.path(
  args$work_dir,
  "savedObjects",
  fds_name,
  "fds-object.RDS"
)

if (file.exists(fds_save_path)) {
    message("Successfully initialized FDS skeleton at: ", fds_save_path)
    # Output path for Hail logs
    cat(fds_save_path, "\n")
} else {
    stop("FDS object was not saved to: ", fds_save_path)
}
