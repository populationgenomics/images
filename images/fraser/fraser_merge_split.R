#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(BiocParallel)
library(summarizedRT) # Useful for rowMaxs/assay handling

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
file.copy(args$fds_path, file.path(save_dir, "fds-object.RDS"))

# 2. Configure Backend
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)
bp <- MulticoreParam(workers = args$nthreads)
register(bp)

# 3. Load Dataset
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)

# 4. Merge individual split count RDS files from the cache
# FRASER automatically looks in: {work_dir}/cache/splitCounts/
message("Merging split counts from cache...")
fds <- countRNAData(fds, recount = FALSE, BPPARAM = bp)

# 5. Extract and Annotate Junctions
# This identifies which junctions exist across the whole cohort
split_counts_se <- asSE(fds, type = "j")
split_ranges <- rowRanges(split_counts_se)

# Explicitly annotate splice sites to get donor/acceptor positions
message("Annotating splice sites...")
split_ranges <- FRASER:::annotateSpliceSite(split_ranges)
saveRDS(split_ranges, "g_ranges_split_counts.RDS")

# 6. Filtering for Non-Split Counting (Optimization)
# We only want to count non-split reads for junctions that actually show up
# in our data to save massive amounts of compute time.
message("Filtering ranges for non-split counting...")
minExpressionInOneSample <- 20
# Use assay(..., "rawCountsJ") to get the counts for junctions
raw_counts <- assay(split_counts_se, "rawCountsJ")
max_count <- rowMaxs(raw_counts)
passed <- max_count >= minExpressionInOneSample

filtered_ranges <- split_ranges[passed, ]
saveRDS(filtered_ranges, "g_ranges_non_split_counts.RDS")

# 7. Extract Splice Site Coordinates
# These are the specific base-pair positions the next jobs will query in the BAMs
message("Extracting splice site coordinates...")
splice_site_coords <- FRASER:::extractSpliceSiteCoordinates(filtered_ranges, fds)
saveRDS(splice_site_coords, "splice_site_coords.RDS")

# 8. Save the updated FDS
saveFraserDataSet(fds)
message("Merge split complete.")
