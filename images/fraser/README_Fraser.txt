The tool is called "Fraser" and it is an R package that can be used to analyze RNA-seq data for splicing events.
It is designed to identify and quantify splicing events in a cohort of samples, and can be used to compare
splicing patterns between different groups of samples.

This image and tool was made to run Cohort-Level analyses on Exome data for the RD team.
When using these scripts and tools, please make sure to have the correct version of R and the necessary packages installed.
They were designed to be used in conjunction with the cpg-flow-rdrnaseq and specifically the fraser stage of that workflow and run in this order

1. fraser_init.R:
Initializes the FRASER dataset for a cohort.
Inputs: --cohort_id, --sample_map (CSV with sample_id and bam_path), --work_dir, --nthreads
What it does: Reads the sample map, creates a FRASER dataset skeleton, and saves it for downstream steps.

2. fraser_count_split.R:
Counts split (spliced) reads for a single sample.
Inputs: --fds_path, --bam_path, --cohort_id, --sample_id, --work_dir, --nthreads
What it does: Loads the dataset, subsets to the sample, counts split reads, and saves results to the cache.

3. fraser_merge_split.R:
Merges split read counts from all samples.
Inputs: --fds_path, --cohort_id, --work_dir, --nthreads
What it does: Loads and merges all split count files, annotates splice sites, and saves the merged data for the cohort.

4. fraser_count_non_split.R:
Counts non-split (unspliced) reads for a single sample.
Inputs: --fds_path, --bam_path, --cohort_id, --sample_id, --coords_path, --work_dir, --nthreads
What it does: Counts non-split reads at splice site coordinates and saves results to the cache.

5. fraser_merge_non_split.R:
Merges non-split read counts and calculates PSI values.
Inputs: --fds_path, --cohort_id, --filtered_ranges_path, --work_dir, --nthreads
What it does: Merges non-split counts for all samples, organizes HDF5 files, and prepares data for analysis.

6. fraser_join_counts.R:
Joins split and non-split counts into a single dataset.
Inputs: --fds_path, --cohort_id, --work_dir, --nthreads
What it does: Checks data integrity, joins split and non-split counts, and ensures all metadata is present for analysis.

7. fraser_analysis.R:
Performs statistical analysis and generates QC plots.
Inputs: --fds_dir, --cohort_id, --pval_cutoff, --delta_psi_cutoff, --min_count, --nthreads
What it does: Loads the final dataset, calculates splicing statistics, applies filters, and generates quality control plots and results.
