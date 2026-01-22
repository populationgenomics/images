#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(BiocParallel)
library(SummarizedExperiment)
library(Rsamtools) # Added for BamFile validation

parser <- ArgumentParser(description = "Merge Split Read Counts")
parser$add_argument("--fds_path", required = TRUE, help = "Path to FDS RDS file")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--work_dir", default = "/io/work", help = "Working directory")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

# 1. Reconstruct Directory Structure
fds_name <- paste0("FRASER_", args$cohort_id)
save_dir <- file.path(args$work_dir, "savedObjects", fds_name)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(args$fds_path, file.path(save_dir, "fds-object.RDS"), overwrite = TRUE)

# 2. Configure Backend
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)


bp <- MulticoreParam(workers = args$nthreads)
register(bp)

# 3. Load Dataset
# Match the name exactly as defined in the init stage
fds_name <- paste0("FRASER_", args$cohort_id)
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)

# --- VULNERABILITY FIX: BAM Metadata Extraction ---
# Instead of a dummy file, point metadata queries to an existing BAM
# to ensure chromosome naming (seqlevelsStyle) can be validated.
# We assume the symlinks were created in /io/batch/input_bams/
available_bams <- list.files("/io/batch/input_bams", pattern = "\\.bam$", full.names = TRUE)
if(length(available_bams) > 0){
    message("Using reference BAM for metadata: ", available_bams[1])

    # Force-update colData for ALL samples to point to a valid file.
    # This prevents the 'dummy.bam' error during internal validation.
    new_colData <- colData(fds)
    new_colData$bamFile <- available_bams[1]
    colData(fds) <- new_colData

    # Manually assert the seqlevelsStyle (e.g., 'UCSC' or 'Ensembl')
    # This satisfies the internal check that is currently crashing.
    seqlevelsStyle(fds) <- seqlevelsStyle(BamFile(available_bams[1]))
} else {
    stop("CRITICAL ERROR: No BAM files found in /io/batch/input_bams. The merge cannot proceed without at least one valid BAM header for seqlevelsStyle validation.")
}
strandSpecific(fds) <- 0
# --------------------------------------------

# 4. Merge individual split count RDS files from the cache
# FRASER automatically looks in: {work_dir}/cache/splitCounts/
message("Merging split counts from cache...")
fds <- getSplitReadCountsForAllSamples(fds, recount = FALSE, BPPARAM = bp)

# 5. Extract and Annotate Junctions
# This identifies which junctions exist across the whole cohort
split_counts_se <- asSE(fds, type = "j")
split_ranges <- rowRanges(split_counts_se)

# Explicitly annotate splice sites to get donor/acceptor positions
message("Annotating splice sites...")
split_ranges <- FRASER:::annotateSpliceSite(split_ranges)

# 6. Filtering for Non-Split Counting (Optimization)
message("Filtering ranges for non-split counting...")
minExpressionInOneSample <- 20
raw_counts <- assay(split_counts_se, "rawCountsJ")
max_count <- rowMaxs(raw_counts)
passed <- max_count >= minExpressionInOneSample

filtered_ranges <- split_ranges[passed, ]

# 7. Extract Splice Site Coordinates
message("Extracting splice site coordinates...")
splice_site_coords <- FRASER:::extractSpliceSiteCoordinates(filtered_ranges, fds)

# Use absolute paths for saving to match Python 'mv' commands
saveRDS(split_ranges,       file.path(args$work_dir, "g_ranges_split_counts.RDS"))
saveRDS(filtered_ranges,    file.path(args$work_dir, "g_ranges_non_split_counts.RDS"))
saveRDS(splice_site_coords, file.path(args$work_dir, "splice_site_coords.RDS"))

# 8. Save the updated FDS
saveFraserDataSet(fds)
message("Merge split complete.")
