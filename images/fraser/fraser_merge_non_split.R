#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(BiocParallel)

parser <- ArgumentParser(description = "Merge Non-Split Counts and Calculate PSI")
parser$add_argument("--fds_path", required = TRUE, help = "Path to localized FDS RDS file")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--filtered_ranges_path", required = TRUE, help = "Path to g_ranges_non_split_counts.RDS")
parser$add_argument("--work_dir", default = "/io/work", help = "Working directory")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

# 1. Reconstruct Directory Structure
fds_name <- paste0("FRASER_", args$cohort_id)
# Where the pipeline put them
old_dir <- file.path(args$work_dir, "savedObjects", "Data_Analysis", "nonSplitCounts")
# Where FRASER expects them
new_dir <- file.path(args$work_dir, "savedObjects", fds_name, "nonSplitCounts")

# --- 2. Move Files to Correct Location ---
if (dir.exists(old_dir)) {
    message("Moving .h5 files from ", old_dir, " to ", new_dir)
    dir.create(new_dir, recursive = TRUE, showWarnings = FALSE)

    h5_files <- list.files(old_dir, pattern = "\\.h5$", full.names = TRUE)
    for (f in h5_files) {
        dest <- file.path(new_dir, basename(f))
        file.rename(f, dest)
    }
} else {
    message("Notice: Source directory ", old_dir, " not found. Checking if files are already in place.")
}

# --- 3. Prepare FDS for Loading ---
# FRASER's loadFraserDataSet looks for the .RDS file in a specific subfolder
save_dir <- file.path(args$work_dir, "savedObjects", fds_name)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(args$fds_path, file.path(save_dir, "fds-object.RDS"), overwrite = TRUE)

# --- 4. Load Dataset ---
# dir is the ROOT workdir, name is the cohort folder name
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)


available_bams <- list.files("/io/batch/input_bams", pattern = "\\.bam$", full.names = TRUE)
if(length(available_bams) > 0){
    ref_bam <- available_bams[1]
    colData(fds)$bamFile <- ref_bam
    seqlevelsStyle(fds) <- seqlevelsStyle(Rsamtools::BamFile(ref_bam))
}
strandSpecific(fds) <- 0

# --- 6. Merge Non-Split Counts ---
message("Merging counts using HDF5 files in: ", new_dir)
non_split_count_ranges <- readRDS(args$filtered_ranges_path)
non_split_counts <- getNonSplitReadCountsForAllSamples(
  fds = fds,
  splitCountRanges = non_split_count_ranges,
  minAnchor = 5,
  recount = FALSE # Crucial: FALSE ensures it uses the .h5 files from the cache
)

# 7. Final Save
# This creates the complete 'fitted' dataset that the analysis script will use

saveRDS(spliceSiteCoords, file.path(args$work_dir, "splice_site_coords.RDS"))

