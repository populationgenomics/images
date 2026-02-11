#!/usr/bin/env Rscript

library(argparse)
library(FRASER)
library(BiocParallel)
library(SummarizedExperiment)
library(HDF5Array)
library(rhdf5)

parser <- ArgumentParser(description = "Merge Non-Split Counts and Calculate PSI")
parser$add_argument("--fds_path", required = TRUE, help = "Path to localized FDS RDS file")
parser$add_argument("--cohort_id", required = TRUE, help = "Cohort ID")
parser$add_argument("--filtered_ranges_path", required = TRUE, help = "Path to splice_site_coords.RDS")
parser$add_argument("--work_dir", default = "/io/work", help = "Working directory")
parser$add_argument("--nthreads", type = "integer", default = 1, help = "Number of threads")
args <- parser$parse_args()

# Configure FRASER to use HDF5 backend
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)

# 1. Setup directories
message("\n[1/6] Setting up directories...")
fds_name <- paste0("FRASER_", args$cohort_id)
save_dir <- file.path(args$work_dir, "savedObjects", fds_name)
out_dir <- file.path(save_dir, "nonSplitCounts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Copy FDS file only if source and destination are different
fds_dest <- file.path(save_dir, "fds-object.RDS")
if(normalizePath(args$fds_path, mustWork = FALSE) != normalizePath(fds_dest, mustWork = FALSE)) {
    file.copy(args$fds_path, fds_dest, overwrite = TRUE)
    message("Copied FDS from: ", args$fds_path)
} else {
    message("FDS already in correct location: ", fds_dest)
}

# 2. Load FDS to get sample information
message("\n[2/6] Loading FraserDataSet...")
fds <- loadFraserDataSet(dir = args$work_dir, name = fds_name)
sample_ids <- samples(fds)
n_samples <- length(sample_ids)
message("Samples: ", n_samples)

# 3. Load the genomic ranges
message("\n[3/6] Loading genomic ranges...")
non_split_ranges <- readRDS(args$filtered_ranges_path)
n_sites <- length(non_split_ranges)
message("Non-split sites: ", n_sites)

# 4. Copy H5 files from cache to the proper output directory
message("\n[4/6] Copying and organizing H5 files...")
cache_dir <- file.path(args$work_dir, "cache", "nonSplicedCounts")
h5_files <- list.files(cache_dir, pattern = "\\.h5$", full.names = TRUE)
message("Found ", length(h5_files), " HDF5 files in cache")

if(length(h5_files) == 0) {
    stop("ERROR: No H5 files found in cache directory: ", cache_dir)
}

# Copy with standardized naming convention
for(i in seq_along(sample_ids)) {
    sid <- sample_ids[i]
    src_file <- h5_files[grepl(sid, h5_files)][1]

    if(is.na(src_file)) {
        warning("No H5 file found for sample: ", sid)
        next
    }

    dest <- file.path(out_dir, paste0("nonSplicedCounts-", sid, ".h5"))
    file.copy(src_file, dest, overwrite = TRUE)
    message("  Copied: ", basename(src_file), " -> ", basename(dest))
}

# 5. Validate H5 files and build combined count matrix
message("\n[5/6] Building combined count matrix (memory-efficient)...")

# Read first file to validate dimensions
first_h5 <- file.path(out_dir, paste0("nonSplicedCounts-", sample_ids[1], ".h5"))
if(!file.exists(first_h5)) {
    stop("ERROR: First H5 file not found: ", first_h5)
}

h5_info <- h5ls(first_h5)
dataset_row <- h5_info[h5_info$name == "nonSplicedCounts", ]
actual_rows <- as.numeric(strsplit(dataset_row$dim, " x ")[[1]][1])

if(actual_rows != n_sites) {
    stop("ERROR: H5 files have ", actual_rows, " rows but ranges have ", n_sites,
         " sites. Mismatch!")
}

message("Validated: H5 dimensions match genomic ranges (", actual_rows, " sites)")

# Create a combined count matrix by reading each sample one at a time
count_matrix <- matrix(0, nrow = n_sites, ncol = n_samples)
colnames(count_matrix) <- sample_ids

total_counts <- 0
for(i in seq_along(sample_ids)) {
    sid <- sample_ids[i]
    h5_file <- file.path(out_dir, paste0("nonSplicedCounts-", sid, ".h5"))

    if(!file.exists(h5_file)) {
        warning("H5 file not found for sample ", sid, ", using zeros")
        next
    }

    message("  Reading [", i, "/", n_samples, "]: ", sid)

    # Read the counts data for this sample only
    counts <- h5read(h5_file, "nonSplicedCounts")
    if(is.matrix(counts)) {
        counts <- counts[, 1]
    }

    count_matrix[, i] <- as.vector(counts)
    total_counts <- total_counts + sum(counts)
}

message("Combined matrix dimensions: ", nrow(count_matrix), " x ", ncol(count_matrix))
message("Total counts: ", total_counts)

if(total_counts == 0) {
    stop("ERROR: The merged non-split counts matrix is empty (all zeros).")
}

# 6. Create HDF5-backed SummarizedExperiment
message("\n[6/6] Creating HDF5-backed SummarizedExperiment...")

# Create SE with in-memory matrix
nonSplitCounts_se <- SummarizedExperiment(
    assays = list(rawCountsSS = count_matrix),
    rowRanges = non_split_ranges,
    colData = DataFrame(sampleID = sample_ids, row.names = sample_ids)
)

message("Non-split SE created with dimensions: ", nrow(nonSplitCounts_se),
        " x ", ncol(nonSplitCounts_se))

# Annotate non-split counts with splice site IDs
message("Annotating non-split counts with splice site IDs...")
nonSplitCounts_se <- FRASER:::annotateSpliceSite(nonSplitCounts_se)

# Save the SummarizedExperiment - this will automatically create HDF5 backend
message("Saving non-split HDF5 SummarizedExperiment...")
saveHDF5SummarizedExperiment(nonSplitCounts_se, dir = out_dir, replace = TRUE)

# Reload it to get the proper HDF5-backed version
message("Reloading HDF5-backed version...")
non_split_counts <- loadHDF5SummarizedExperiment(dir = out_dir)

# Final verification
if(!"rawCountsSS" %in% assayNames(non_split_counts)) {
    stop("ERROR: 'rawCountsSS' assay not found in reloaded object.")
}

message("Verification: Reloaded object has ", sum(assay(non_split_counts, "rawCountsSS")), " total counts")

# 7. Configure Backend for FDS
message("\nConfiguring FRASER backend...")
bp <- MulticoreParam(workers = args$nthreads)
register(bp)

# Set BAM file metadata (for validation purposes)
available_bams <- list.files("/io/batch/input_bams", pattern = "\\.bam$", full.names = TRUE)
if(length(available_bams) > 0) {
    ref_bam <- available_bams[1]
    message("Validating metadata against reference BAM: ", ref_bam)
    colData(fds)$bamFile <- ref_bam
    seqlevelsStyle(fds) <- seqlevelsStyle(Rsamtools::BamFile(ref_bam))
} else {
    warning("WARNING: No BAM files found in /io/batch/input_bams.")
}
strandSpecific(fds) <- 2

# 8. Update FDS object with the merged non-split counts
message("\nUpdating FDS object with merged non-split counts...")
nonSplicedReads(fds) <- non_split_counts
fds <- saveFraserDataSet(fds)

message("\n=== Non-split merge complete! ===")
message("Final stats:")
message("  - Sites: ", nrow(non_split_counts))
message("  - Samples: ", ncol(non_split_counts))
message("  - Total counts: ", sum(assay(non_split_counts, "rawCountsSS")))
