#!/usr/bin/env Rscript
library(argparse)
library(FRASER)
library(SummarizedExperiment)
library(HDF5Array)

parser <- ArgumentParser(description = "Join Split and Non-Split Counts")
parser$add_argument("--fds_path", required = TRUE)
parser$add_argument("--cohort_id", required = TRUE)
parser$add_argument("--work_dir", default = "/io/work")
parser$add_argument("--nthreads", type = "integer", default = 1)
args <- parser$parse_args()

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
    # Are all junction startIDs present in the Splice Site map?
    missing_starts <- !all(mcols(fds, type="j")$startID %in% mcols(fds, type="ss")$spliceSiteID)
    if(missing_starts) {
        stop("CRITICAL ERROR: Junction startIDs do not match SpliceSiteIDs. Map is broken.")
    }
    message("SUCCESS: Junction IDs are correctly mapped to Splice Sites.")

    # 3. Check HDF5 Backend connectivity
    # Ensure the pointers to the .h5 files are valid and readable
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
message("Loading Fraser Data Set from RDS...")
fds <- readRDS(args$fds_path)
workingDir(fds) <- saveDir

# 2. Load split counts
message("Loading split counts...")
split_se_dir <- file.path(saveDir, "splitCounts")
if(!dir.exists(split_se_dir)){
    stop(paste("Missing splitCounts directory at:", split_se_dir))
}
splitCounts_se <- loadHDF5SummarizedExperiment(dir = split_se_dir)

# 3. Load merged non-split counts
message("Loading merged non-split counts...")
merged_non_split_dir <- file.path(saveDir, "nonSplitCounts")
if(!dir.exists(merged_non_split_dir)){
    stop(paste("Missing merged non-split counts directory at:", merged_non_split_dir))
}
nonSplitCounts_se <- loadHDF5SummarizedExperiment(dir = merged_non_split_dir)

# 4. Annotate the Split Counts
# This is the missing link. It generates the spliceSiteID mapping
# that calculatePSIValues needs later.
message("Annotating splice sites...")
splitCounts_se <- FRASER:::annotateSpliceSite(splitCounts_se)

# 5. Add counts to FRASER object
message("Joining assays into FDS object...")
fds <- addCountsToFraserDataSet(
  fds = fds,
  splitCounts = splitCounts_se,
  nonSplitCounts = nonSplitCounts_se
)

# --- Call the check before finishing ---
check_fds_integrity(fds)

# This will populate the 'ss' (splice site) metadata correctly
fds <- calculatePSIValues(fds, types="jaccard")

# --- Call the check before finishing ---
check_fds_integrity(fds)
# 5. Save final FRASER object
message("Saving final integrated FDS...")
fds <- saveFraserDataSet(fds)

message("FRASER join complete.")
