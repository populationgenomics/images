#!/usr/bin/env Rscript
library(argparse)
library(FRASER)
library(SummarizedExperiment)
library(HDF5Array)
library(BiocParallel)

parser <- ArgumentParser(description = "Join Split and Non-Split Counts")
parser$add_argument("--fds_path", required = TRUE)
parser$add_argument("--cohort_id", required = TRUE)
parser$add_argument("--work_dir", default = "/io/work")
parser$add_argument("--nthreads", type = "integer", default = 1)
args <- parser$parse_args()

# Configure FRASER to use HDF5 backend
options("FRASER.maxSamplesNoHDF5" = 0)
options("FRASER.maxJunctionsNoHDF5" = -1)

# Setup parallel processing
bp <- MulticoreParam(workers = args$nthreads)
register(bp)

check_fds_integrity <- function(fds) {
    message("\n--- Running FDS Integrity Check ---")

    # 1. Check for basic ID existence
    j_cols <- colnames(mcols(fds, type="j"))
    ss_cols <- colnames(mcols(fds, type="ss"))

    has_ids <- all(c("startID", "endID") %in% j_cols) && ("spliceSiteID" %in% ss_cols)

    if(has_ids) {
        message("SUCCESS: SpliceSiteIDs found in both Junctions and SpliceSites.")
    } else {
        stop("CRITICAL ERROR: SpliceSiteID metadata is missing. Analysis will fail.")
    }

    # 2. Check for ID alignment (The "by.y" error test)
    missing_starts <- !all(mcols(fds, type="j")$startID %in% mcols(fds, type="ss")$spliceSiteID)
    if(missing_starts) {
        stop("CRITICAL ERROR: Junction startIDs do not match SpliceSiteIDs. Map is broken.")
    }
    message("SUCCESS: Junction IDs are correctly mapped to Splice Sites.")

    # 3. Check HDF5 Backend connectivity
    tryCatch({
        test_val <- as.matrix(counts(fds, type="j")[1, 1])
        message("SUCCESS: HDF5 backends are reachable and readable.")
    }, error = function(e) {
        stop("CRITICAL ERROR: HDF5 files are missing or paths are broken: ", e$message)
    })

    # 4. Check for Non-Split Counts
    if("rawCountsSS" %in% assayNames(fds)) {
        message("SUCCESS: Non-split counts (rawCountsSS) are present.")
    } else {
        warning("WARNING: Non-split counts are missing. Jaccard/PSI calculation may fail.")
    }

    message("--- Integrity Check Passed ---\n")
}

fds_name <- paste0("FRASER_", args$cohort_id)
saveDir <- file.path(args$work_dir, "savedObjects", fds_name)

# 1. Load the FDS object
message("\n[1/7] Loading Fraser Data Set from RDS...")
fds <- readRDS(args$fds_path)
workingDir(fds) <- saveDir

# 2. Load the master splice site coordinates that were used to create non-split counts
message("\n[2/7] Loading master splice site coordinates...")
splice_coords_path <- file.path(dirname(dirname(saveDir)), "splice_site_coords.RDS")
if(!file.exists(splice_coords_path)) {
    stop("ERROR: splice_site_coords.RDS not found at: ", splice_coords_path)
}
master_splice_coords <- readRDS(splice_coords_path)
message("Master coordinates loaded: ", length(master_splice_coords), " sites")

# 3. Load split counts
message("\n[3/7] Loading split counts...")
split_se_dir <- file.path(saveDir, "splitCounts")
if(!dir.exists(split_se_dir)){
    stop(paste("Missing splitCounts directory at:", split_se_dir))
}
splitCounts_se <- loadHDF5SummarizedExperiment(dir = split_se_dir)
message("Split counts dimensions: ", nrow(splitCounts_se), " junctions x ", ncol(splitCounts_se), " samples")

# Annotate split counts with splice site IDs
# NOTE: annotateSpliceSite() automatically extracts splice sites from junctions
# and generates consistent startID/endID mappings
message("Annotating split counts with splice site IDs...")
splitCounts_se <- FRASER:::annotateSpliceSite(splitCounts_se)

# 4. Load merged non-split counts
message("\n[4/7] Loading merged non-split counts...")
merged_non_split_dir <- file.path(saveDir, "nonSplitCounts")
if(!dir.exists(merged_non_split_dir)){
    stop(paste("Missing merged non-split counts directory at:", merged_non_split_dir))
}
nonSplitCounts_se <- loadHDF5SummarizedExperiment(dir = merged_non_split_dir)
message("Non-split counts dimensions: ", nrow(nonSplitCounts_se), " sites x ", ncol(nonSplitCounts_se), " samples")

# CRITICAL FIX: Always re-annotate after loading to ensure consistency
# The saved HDF5 object might have stale or incompatible IDs
message("Re-annotating non-split counts with splice site IDs...")
nonSplitCounts_se <- FRASER:::annotateSpliceSite(nonSplitCounts_se)

# 5. Verify ID consistency before joining
message("\n[5/7] Verifying ID consistency...")
split_start_ids <- na.omit(unique(mcols(splitCounts_se)$startID))
split_end_ids <- na.omit(unique(mcols(splitCounts_se)$endID))
nonsplit_ids <- na.omit(unique(mcols(nonSplitCounts_se)$spliceSiteID))

message("  Split startIDs: ", length(split_start_ids))
message("  Split endIDs: ", length(split_end_ids))
message("  Non-split spliceSiteIDs: ", length(nonsplit_ids))

# Check overlap - both split start/end should be subsets of nonsplit
missing_start <- sum(!split_start_ids %in% nonsplit_ids)
missing_end <- sum(!split_end_ids %in% nonsplit_ids)

if(missing_start > 0 || missing_end > 0) {
    message("NOTE: Some split junction IDs are not in non-split IDs:")
    message("  Missing startIDs: ", missing_start, " out of ", length(split_start_ids))
    message("  Missing endIDs: ", missing_end, " out of ", length(split_end_ids))
    message("  This is expected if junctions reference splice sites not in the filtered set")
}

# 6. Add counts to FRASER object
message("\n[6/7] Joining split and non-split counts into FDS object...")
fds <- addCountsToFraserDataSet(
  fds = fds,
  splitCounts = splitCounts_se,
  nonSplitCounts = nonSplitCounts_se
)

message("Counts successfully joined!")
message("  - Split junctions: ", nrow(counts(fds, type = "j")))
message("  - Splice sites: ", nrow(counts(fds, type = "ss")))

# Run integrity check
check_fds_integrity(fds)

# 7. Calculate PSI values
message("\n[7/7] Calculating PSI values...")
fds <- calculatePSIValues(fds, types = c("psi3", "psi5", "jaccard"))

message("PSI values calculated successfully!")
message("Available assays: ", paste(assayNames(fds), collapse = ", "))

# Final integrity check
check_fds_integrity(fds)

# 8. Save final FRASER object
message("\nSaving final integrated FDS...")
fds <- saveFraserDataSet(fds)

message("\n=== FRASER join complete! ===")
message("FDS object saved in: ", workingDir(fds))
message("\nNext steps for analysis:")
message("  1. Filter junctions: fds <- filterExpressionAndVariability(fds)")
message("  2. Fit model: fds <- FRASER(fds)")
message("  3. Get results: results(fds)")
