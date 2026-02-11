#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(BiocParallel)
library(SummarizedExperiment)
library(Rsamtools)
library(DelayedMatrixStats)
library(HDF5Array)


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
    stop("CRITICAL ERROR: No BAM files found in /io/batch/input_bams.")
}
strandSpecific(fds) <- 2
# --------------------------------------------
minExpressionInOneSample <- 20
# 4. Merge individual split count RDS files from the cache
# FRASER automatically looks in: {work_dir}/cache/splitCounts/
message("Merging split counts from cache...")
splitCounts <- getSplitReadCountsForAllSamples(fds=fds, recount=FALSE)

# --- CRITICAL FIX: Annotate the object ITSELF, not a temporary variable ---
# This generates startID and endID and attaches them to the rowRanges
message("Generating Splice Site IDs...")
rowRanges(splitCounts) <- FRASER:::annotateSpliceSite(rowRanges(splitCounts))

# 5. Now save the Annotated SE
message("Saving merged split counts SummarizedExperiment with IDs...")
split_counts_dir <- file.path(save_dir, "splitCounts")
dir.create(split_counts_dir, recursive = TRUE, showWarnings = FALSE)
saveHDF5SummarizedExperiment(splitCounts, dir = split_counts_dir, replace = TRUE)

# 6. Extract coordinates for the next pipeline steps
# We use the rowRanges that NOW HAVE the startID/endID columns
splitCountRangesNonFilt <- rowRanges(splitCounts)

maxCount <- rowMaxs(assay(splitCounts, "rawCountsJ"))
passed <- maxCount >= minExpressionInOneSample
# extract granges after filtering
splitCountRanges <- splitCountRangesNonFilt[passed,]
spliceSiteCoords <- FRASER:::extractSpliceSiteCoordinates(splitCountRanges)


# Use absolute paths for saving to match Python 'mv' commands
#This is slightly confusing, but the filtered granges will be used to annotate non_split counts
saveRDS(splitCountRangesNonFilt,       file.path(args$work_dir, "g_ranges_split_counts.RDS"))
saveRDS(splitCountRanges,    file.path(args$work_dir, "g_ranges_non_split_counts.RDS"))
saveRDS(spliceSiteCoords, file.path(args$work_dir, "splice_site_coords.RDS"))
