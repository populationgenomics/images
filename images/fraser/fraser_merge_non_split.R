#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(BiocParallel)
library(SummarizedExperiment)
library(HDF5Array)

parser <- ArgumentParser(description = "Merge Non-Split Counts and Calculate PSI")
parser$add_argument("--fds_path", required = TRUE, help = "Path to localized FDS RDS file")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--filtered_ranges_path", required = TRUE, help = "Path to g_ranges_non_split_counts.RDS")
parser$add_argument("--work_dir", default = "/io/work", help = "Working directory")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

# 1. Load the ranges first so they are available for all steps
non_split_count_ranges <- readRDS(args$filtered_ranges_path)
fds_name <- paste0("FRASER_", args$cohort_id)
save_dir <- file.path(args$work_dir, "savedObjects", fds_name)
out_dir <- file.path(save_dir, "nonSplitCounts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(args$fds_path, file.path(save_dir, "fds-object.RDS"), overwrite = TRUE)

# 2. Load FDS
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)

# 3. Copy H5 files from cache to the proper output directory
cache_dir <- file.path(args$work_dir, "cache", "nonSplicedCounts")
h5_files <- list.files(cache_dir, pattern = "\\.h5$", full.names = TRUE)
message("Copying ", length(h5_files), " HDF5 files from cache to output directory...")
for(f in h5_files) {
    file.copy(f, file.path(out_dir, basename(f)), overwrite = TRUE)
}

anchor_se <- SummarizedExperiment(
    assays = list(rawCountsSS = matrix(0,
                                       nrow = length(non_split_count_ranges),
                                       ncol = length(samples(fds)),
                                       dimnames = list(NULL, samples(fds)))),
    rowRanges = non_split_count_ranges
)
# This creates the se.rds that FRASER's loadHDF5SummarizedExperiment needs
saveHDF5SummarizedExperiment(anchor_se, dir = out_dir, replace = TRUE, prefix = "rawCountsSS")

# 4. Configure Backend
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)
bp <- MulticoreParam(workers = args$nthreads)
register(bp)


available_bams <- list.files("/io/batch/input_bams", pattern = "\\.bam$", full.names = TRUE)
if(length(available_bams) > 0){
    ref_bam <- available_bams[1]
    message("Validating metadata against reference BAM: ", ref_bam)
    colData(fds)$bamFile <- ref_bam
    seqlevelsStyle(fds) <- seqlevelsStyle(Rsamtools::BamFile(ref_bam))
} else {
    stop("CRITICAL ERROR: No BAM files found in /io/batch/input_bams.")
}
strandSpecific(fds) <- 2

# --- 6. Merge Non-Split Counts ---
message("Merging counts using HDF5 files in: ", out_dir)
non_split_counts <- getNonSplitReadCountsForAllSamples(
  fds = fds,
  splitCountRanges = non_split_count_ranges,
  minAnchor = 5,
  recount = FALSE # Crucial: FALSE ensures it uses the .h5 files from the cache
)

# CRITICAL ADDITION: Annotate the non-split counts object!
# This assigns 'spliceSiteID' to the rowRanges of your non-split object.
message("Annotating non-split splice site IDs...")
non_split_counts <- FRASER:::annotateSpliceSite(non_split_counts)

# --- Verification Block ---
# 1. Check if the assay exists
if(!"rawCountsSS" %in% assayNames(non_split_counts)){
    stop("Merge failed: 'rawCountsSS' assay not found in merged object.")
}

# 2. Check for zeros (Common sign of a path/naming mismatch)
total_counts <- sum(assay(non_split_counts, "rawCountsSS"))
if(total_counts == 0){
    stop("Merge Error: The merged non-split counts matrix is empty (all zeros).")
}

message("Merge verified: Found ", total_counts, " total site coverage counts.")

# 8. Save the merged non-split counts to the correct location for join step
message("Saving merged non-split counts to: ", out_dir)
saveHDF5SummarizedExperiment(
    non_split_counts,
    dir = out_dir,
    replace = TRUE
)

# 9. Update FDS object and save
nonSplicedReads(fds) <- non_split_counts
fds <- saveFraserDataSet(fds)

message("Non-split merge complete!")
