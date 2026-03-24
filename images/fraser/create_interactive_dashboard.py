#!/usr/bin/env python3
"""
Create an interactive dashboard with volcano plot, manhattan plot, and IGV browser.

Supports FRASER and OUTRIDER results with optional TABIX indexing for efficient
region-based queries. Optionally annotates FRASER results with proximity to rare
variants from a Hail table.

Usage:
    python create_interactive_dashboard.py --fraser fraser_results.csv --output dashboard.html
    python create_interactive_dashboard.py --fraser fraser.tsv.gz --outrider outrider.tsv.gz --output dashboard.html
    python create_interactive_dashboard.py --fraser fraser.csv --variant-table variants.ht --output dashboard.html

Input file formats:
    - CSV files: Will be converted to tabix-indexed format automatically
    - .gz files with .tbi index: Used directly for tabix queries
    - Hail tables (.ht): Used for variant proximity annotation

Required FRASER columns: seqnames, start, end, pValue, padjust, deltaPsi, psiValue, type, hgncSymbol
Required OUTRIDER columns: seqnames, start, end, pValue, padjust, (optional: zScore, log2fc, geneID/hgncSymbol)

Coordinate system: 1-based, inclusive (BED-like but 1-based for consistency with R/bioconductor)
"""

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader

# Try to import pysam for tabix support
try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False
    print("Warning: pysam not installed. Tabix queries will use subprocess fallback.")

# Try to import hail for variant table support
try:
    import hail as hl
    HAS_HAIL = True
except ImportError:
    HAS_HAIL = False
    print("Warning: hail not installed. Variant proximity annotations will not be available.")


# =============================================================================
# Constants and Configuration
# =============================================================================

FRASER_REQUIRED_COLUMNS = ['seqnames', 'start', 'end', 'pValue', 'padjust',
                           'deltaPsi', 'psiValue', 'type', 'sampleID']
FRASER_OPTIONAL_COLUMNS = ['hgncSymbol', 'counts', 'totalCounts', 'nonsplitCounts']

OUTRIDER_REQUIRED_COLUMNS = ['pValue', 'padjust', 'sampleID']
OUTRIDER_OPTIONAL_COLUMNS = ['zScore', 'l2fc', 'log2fc', 'geneID', 'hgncSymbol', 'seqnames', 'start', 'end']

# Chromosome sort order for Manhattan plot
CHR_ORDER = {
    **{str(i): i for i in range(1, 23)},
    **{f'chr{i}': i for i in range(1, 23)},
    'X': 23, 'chrX': 23,
    'Y': 24, 'chrY': 24,
    'M': 25, 'MT': 25, 'chrM': 25, 'chrMT': 25
}


# =============================================================================
# CPG to Family ID Mapping
# =============================================================================

def load_cpg_to_family_mapping(mapping_file: str) -> dict:
    """
    Load CPG ID to Family ID mapping from project summary CSV.

    Args:
        mapping_file: Path to rdnow-export-project-summary CSV file

    Returns:
        Dictionary mapping CPG IDs (sequencing_group.id) to family.external_ids
    """
    print(f"Loading CPG to Family mapping from {mapping_file}...")

    df = pd.read_csv(mapping_file)

    # Create mapping from sequencing_group.id to family.external_ids
    mapping = {}
    for _, row in df.iterrows():
        cpg_id = row['sequencing_group.id']
        family_id = row['family.external_ids']


def load_bam_mapping(mapping_file: str) -> dict:
    """
    Load sample ID to BAM file mapping from TSV.

    Args:
        mapping_file: Path to TSV file with columns: sampleID, bam_url, bai_url

    Returns:
        Dictionary mapping sample IDs to BAM info dicts
    """
    print(f"Loading BAM mapping from {mapping_file}...")

    df = pd.read_csv(mapping_file, sep='\t')

    # Validate required columns
    required_cols = ['sampleID', 'bam_url']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"BAM mapping file must contain '{col}' column")

    # Create mapping
    mapping = {}
    for _, row in df.iterrows():
        sample_id = str(row['sampleID'])
        bam_url = str(row['bam_url'])
        bai_url = str(row.get('bai_url', bam_url + '.bai'))

        mapping[sample_id] = {
            'bam_url': bam_url,
            'bai_url': bai_url
        }

    print(f"  Loaded BAM paths for {len(mapping)} samples")
    return mapping


def load_cpg_to_family_mapping(mapping_file: str) -> dict:
    """
    Load CPG ID to Family ID mapping from project summary CSV.

    Args:
        mapping_file: Path to rdnow-export-project-summary CSV file

    Returns:
        Dictionary mapping CPG IDs (sequencing_group.id) to family.external_ids
    """
    print(f"Loading CPG to Family mapping from {mapping_file}...")

    df = pd.read_csv(mapping_file)

    # Create mapping from sequencing_group.id to family.external_ids
    mapping = {}
    for _, row in df.iterrows():
        cpg_id = row['sequencing_group.id']
        family_id = row['family.external_ids']
        mapping[cpg_id] = family_id

    unique_families = len(set(mapping.values()))
    print(f"  Loaded {len(mapping)} CPG IDs mapping to {unique_families} families")

    return mapping


def add_family_ids(df: pd.DataFrame, cpg_to_family: dict) -> pd.DataFrame:
    """
    Add family ID column to dataframe based on sampleID (CPG ID).

    Args:
        df: DataFrame with 'sampleID' column
        cpg_to_family: Mapping dictionary from load_cpg_to_family_mapping

    Returns:
        DataFrame with added 'familyID' column
    """
    df = df.copy()
    df['familyID'] = df['sampleID'].map(cpg_to_family).fillna('Unknown')

    family_counts = df['familyID'].value_counts()
    print(f"  Mapped {len(df)} variants to {len(family_counts)} families")

    return df


# =============================================================================
# Hail Variant Loading and Proximity Annotation
# =============================================================================

def load_hail_variants(hail_table_path: str) -> pd.DataFrame:
    """
    Load variants from a Hail table and extract genomic positions.

    Args:
        hail_table_path: Path to Hail table (.ht)

    Returns:
        DataFrame with columns: chr, pos, ref, alt (and any additional fields)
    """
    if not HAS_HAIL:
        raise RuntimeError("Hail is not installed. Cannot load Hail variant table.")

    print(f"Loading Hail variant table from {hail_table_path}...")

    # Initialize Hail if not already done
    try:
        hl.init(quiet=True, log='/tmp/hail.log')
    except Exception:
        pass  # Already initialized

    # Load the Hail table
    ht = hl.read_table(hail_table_path)

    # Extract locus and alleles information
    # Assuming the table has a 'locus' field and 'alleles' field (standard Hail format)
    if 'locus' not in ht.row:
        raise ValueError("Hail table must have a 'locus' field")

    # Select relevant fields
    ht_selected = ht.select(
        chr=ht.locus.contig,
        pos=ht.locus.position,
        ref=ht.alleles[0] if 'alleles' in ht.row else '',
        alt=ht.alleles[1] if 'alleles' in ht.row else ''
    )

    # Convert to pandas DataFrame
    df = ht_selected.to_pandas()

    print(f"  Loaded {len(df)} variants from Hail table")
    print(f"  Chromosomes: {sorted(df['chr'].unique())}")

    return df


def find_closest_variant(fraser_df: pd.DataFrame, variants_df: pd.DataFrame,
                         max_distance: int = 1000000) -> pd.DataFrame:
    """
    Annotate Fraser results with the closest rare variant of interest.

    Args:
        fraser_df: DataFrame with Fraser results
        variants_df: DataFrame with variants (chr, pos, ref, alt)
        max_distance: Maximum distance to consider (default 1MB)

    Returns:
        DataFrame with additional columns:
        - closest_variant_chr
        - closest_variant_pos
        - closest_variant_ref
        - closest_variant_alt
        - closest_variant_distance
    """
    print("\nAnnotating Fraser results with closest variants...")

    fraser_annotated = fraser_df.copy()

    # Initialize annotation columns
    fraser_annotated['closest_variant_chr'] = None
    fraser_annotated['closest_variant_pos'] = None
    fraser_annotated['closest_variant_ref'] = None
    fraser_annotated['closest_variant_alt'] = None
    fraser_annotated['closest_variant_distance'] = None

    # Normalize chromosome names in both dataframes
    fraser_annotated['chr_normalized'] = fraser_annotated['seqnames'].apply(normalize_chromosome)
    variants_df_copy = variants_df.copy()
    variants_df_copy['chr_normalized'] = variants_df_copy['chr'].apply(normalize_chromosome)

    # Group variants by chromosome for efficient lookup
    variants_by_chr = {}
    for chr_name in variants_df_copy['chr_normalized'].unique():
        chr_variants = variants_df_copy[variants_df_copy['chr_normalized'] == chr_name].copy()
        # Sort by position for efficient searching
        chr_variants = chr_variants.sort_values('pos')
        variants_by_chr[chr_name] = chr_variants

    print(f"  Processing {len(fraser_annotated)} Fraser results...")

    # For each Fraser result, find the closest variant
    for idx, row in fraser_annotated.iterrows():
        chr_name = row['chr_normalized']
        fraser_start = row['start']
        fraser_end = row['end']
        fraser_middle = (fraser_start + fraser_end) / 2

        # Get variants on the same chromosome
        if chr_name not in variants_by_chr:
            continue

        chr_variants = variants_by_chr[chr_name]

        # Find variants within max_distance
        nearby_variants = chr_variants[
            (chr_variants['pos'] >= fraser_start - max_distance) &
            (chr_variants['pos'] <= fraser_end + max_distance)
        ]

        if len(nearby_variants) == 0:
            continue

        # Calculate distances to Fraser region middle point
        distances = np.abs(nearby_variants['pos'] - fraser_middle)
        min_dist_idx = distances.idxmin()
        closest_variant = nearby_variants.loc[min_dist_idx]

        # Annotate
        fraser_annotated.at[idx, 'closest_variant_chr'] = closest_variant['chr']
        fraser_annotated.at[idx, 'closest_variant_pos'] = int(closest_variant['pos'])
        fraser_annotated.at[idx, 'closest_variant_ref'] = closest_variant['ref']
        fraser_annotated.at[idx, 'closest_variant_alt'] = closest_variant['alt']
        fraser_annotated.at[idx, 'closest_variant_distance'] = int(distances.loc[min_dist_idx])

    # Drop temporary column
    fraser_annotated.drop(columns=['chr_normalized'], inplace=True)

    annotated_count = fraser_annotated['closest_variant_pos'].notna().sum()
    print(f"  Annotated {annotated_count} Fraser results with closest variants")
    print(f"  ({annotated_count / len(fraser_annotated) * 100:.1f}% of total)")

    return fraser_annotated


# =============================================================================
# Chromosome Utilities
# =============================================================================

def normalize_chromosome(chrom: str) -> str:
    """Normalize chromosome name to include 'chr' prefix."""
    chrom = str(chrom)
    if not chrom.startswith('chr'):
        return f'chr{chrom}'
    return chrom


def chr_sort_key(chr_name: str) -> int:
    """Get sort key for chromosome ordering."""
    return CHR_ORDER.get(chr_name, CHR_ORDER.get(str(chr_name).replace('chr', ''), 26))


# =============================================================================
# Data Validation
# =============================================================================

def validate_dataframe(df: pd.DataFrame, required_cols: list, optional_cols: list, name: str) -> None:
    """Validate that a DataFrame has required columns."""
    missing = set(required_cols) - set(df.columns)
    if missing:
        raise ValueError(
            f"{name} data is missing required columns: {missing}\n"
            f"Required: {required_cols}\n"
            f"Found: {list(df.columns)}"
        )

    # Check for optional columns that are present
    present_optional = set(optional_cols) & set(df.columns)
    if present_optional:
        print(f"  Found optional columns: {present_optional}")


# =============================================================================
# Common Results Identification
# =============================================================================

def find_common_significant(df_fraser: pd.DataFrame, df_outrider: pd.DataFrame,
                           pval_threshold: float = 0.05, deltapsi_threshold: float = 0.2) -> tuple:
    """Find genes that are significant in both FRASER and OUTRIDER"""
    print("\nFinding common significant results...")

    # FRASER significant: meets p-value and deltaPSI thresholds
    fraser_sig = df_fraser[
        (df_fraser['padjust'] <= pval_threshold) &
        (np.abs(df_fraser['deltaPsi']) >= deltapsi_threshold)
    ].copy()

    # OUTRIDER significant: meets p-value and z-score thresholds (if zScore column exists)
    if 'zScore' in df_outrider.columns:
        outrider_sig = df_outrider[
            (df_outrider['padjust'] <= pval_threshold) &
            (np.abs(df_outrider['zScore']) >= 2)
        ].copy()
    else:
        outrider_sig = df_outrider[df_outrider['padjust'] <= pval_threshold].copy()

    print(f"  FRASER significant: {len(fraser_sig)}")
    print(f"  OUTRIDER significant: {len(outrider_sig)}")

    # Strategy: Match by sample and gene
    common_by_sample = []
    fraser_samples = set(fraser_sig['sampleID'])
    outrider_samples = set(outrider_sig['sampleID'])
    common_samples = fraser_samples & outrider_samples

    print(f"  Samples with both FRASER and OUTRIDER hits: {len(common_samples)}")

    # For each common sample, find overlapping genes
    for sample in common_samples:
        fraser_sample_genes = fraser_sig[fraser_sig['sampleID'] == sample]
        outrider_sample_genes = outrider_sig[outrider_sig['sampleID'] == sample]

        for _, fraser_row in fraser_sample_genes.iterrows():
            for _, outrider_row in outrider_sample_genes.iterrows():
                # Check if gene names overlap
                fraser_gene = str(fraser_row.get('hgncSymbol', '')).upper() if pd.notna(fraser_row.get('hgncSymbol')) else ''
                outrider_gene = str(outrider_row.get('geneID', outrider_row.get('hgncSymbol', ''))).upper()

                # Consider it a match if the HGNC symbol appears in the gene ID or vice versa
                is_match = False
                if fraser_gene and (
                    fraser_gene in outrider_gene or
                    any(fraser_gene in part for part in outrider_gene.split(';'))
                ):
                    is_match = True

                # If we have a match, add to common results
                if is_match:
                    common_entry = {
                        'gene': fraser_row.get('hgncSymbol', 'NA'),
                        'ensembl_id': outrider_row.get('geneID', outrider_row.get('hgncSymbol', 'NA')),
                        'sampleID': sample,
                        'chr': fraser_row.get('seqnames', 'NA'),
                        'start': fraser_row.get('start', 0),
                        'end': fraser_row.get('end', 0),
                        'fraser_pvalue': fraser_row['pValue'],
                        'fraser_padjust': fraser_row['padjust'],
                        'fraser_deltaPsi': fraser_row['deltaPsi'],
                        'fraser_psiValue': fraser_row.get('psiValue', 'NA'),
                        'fraser_type': fraser_row.get('type', 'NA'),
                        'outrider_pvalue': outrider_row['pValue'],
                        'outrider_padjust': outrider_row['padjust'],
                    }

                    # Add optional OUTRIDER columns
                    if 'zScore' in outrider_row:
                        common_entry['outrider_zScore'] = outrider_row['zScore']
                    if 'l2fc' in outrider_row:
                        common_entry['outrider_l2fc'] = outrider_row['l2fc']
                    elif 'log2fc' in outrider_row:
                        common_entry['outrider_l2fc'] = outrider_row['log2fc']

                    # Add metadata if available
                    if 'familyID' in fraser_row.index:
                        common_entry['familyID'] = fraser_row['familyID']

                    common_by_sample.append(common_entry)

    # Also add entries where sample matches but gene doesn't necessarily match
    # (to show samples with both types of aberrations even if different genes)
    for sample in common_samples:
        fraser_sample_genes = fraser_sig[fraser_sig['sampleID'] == sample]
        outrider_sample_genes = outrider_sig[outrider_sig['sampleID'] == sample]

        # Get the most significant from each
        if len(fraser_sample_genes) > 0 and len(outrider_sample_genes) > 0:
            fraser_top = fraser_sample_genes.nsmallest(1, 'pValue').iloc[0]
            outrider_top = outrider_sample_genes.nsmallest(1, 'pValue').iloc[0]

            # Check if this combination already exists
            already_exists = any(
                r['sampleID'] == sample and
                str(r['gene']) == str(fraser_top.get('hgncSymbol')) and
                str(r['ensembl_id']) == str(outrider_top.get('geneID', outrider_top.get('hgncSymbol')))
                for r in common_by_sample
            )

            if not already_exists:
                common_entry = {
                    'gene': fraser_top.get('hgncSymbol', 'NA'),
                    'ensembl_id': outrider_top.get('geneID', outrider_top.get('hgncSymbol', 'NA')),
                    'sampleID': sample,
                    'chr': fraser_top.get('seqnames', 'NA'),
                    'start': fraser_top.get('start', 0),
                    'end': fraser_top.get('end', 0),
                    'fraser_pvalue': fraser_top['pValue'],
                    'fraser_padjust': fraser_top['padjust'],
                    'fraser_deltaPsi': fraser_top['deltaPsi'],
                    'fraser_psiValue': fraser_top.get('psiValue', 'NA'),
                    'fraser_type': fraser_top.get('type', 'NA'),
                    'outrider_pvalue': outrider_top['pValue'],
                    'outrider_padjust': outrider_top['padjust'],
                }

                # Add optional OUTRIDER columns
                if 'zScore' in outrider_top:
                    common_entry['outrider_zScore'] = outrider_top['zScore']
                if 'l2fc' in outrider_top:
                    common_entry['outrider_l2fc'] = outrider_top['l2fc']
                elif 'log2fc' in outrider_top:
                    common_entry['outrider_l2fc'] = outrider_top['log2fc']

                # Add metadata if available
                if 'familyID' in fraser_top.index:
                    common_entry['familyID'] = fraser_top['familyID']

                common_by_sample.append(common_entry)

    df_common = pd.DataFrame(common_by_sample)
    print(f"  Common significant results: {len(df_common)}")
    print(f"  (Samples showing aberrations in both splicing and expression)")

    return df_common, fraser_sig, outrider_sig


def prepare_table_data(df: pd.DataFrame, columns_subset: Optional[list] = None) -> pd.DataFrame:
    """Prepare data for searchable tables with formatted numeric values"""
    if columns_subset:
        available_columns = [col for col in columns_subset if col in df.columns]
    else:
        available_columns = df.columns.tolist()

    table_df = df[available_columns].copy()

    # Format numeric columns
    for col in table_df.columns:
        if 'pValue' in col or 'padjust' in col or 'pvalue' in col:
            table_df[col] = table_df[col].apply(
                lambda x: f"{x:.2e}" if pd.notna(x) and isinstance(x, (int, float)) else str(x)
            )
        elif 'psi' in col.lower() or 'Score' in col or 'score' in col.lower() or 'l2fc' in col or 'deltaPsi' in col:
            table_df[col] = table_df[col].apply(
                lambda x: f"{x:.3f}" if pd.notna(x) and isinstance(x, (int, float)) else str(x)
            )

    # Fill NaN values
    table_df = table_df.fillna('NA')

    return table_df


# =============================================================================
# Tabix Preparation Functions
# =============================================================================

def is_tabix_indexed(filepath: str) -> bool:
    """Check if a file is already tabix-indexed (.gz with .tbi)."""
    if not filepath.endswith('.gz'):
        return False
    tbi_path = filepath + '.tbi'
    return os.path.exists(tbi_path)


def csv_to_tabix(csv_path: str, output_dir: Optional[str] = None) -> str:
    """
    Convert a CSV file to a sorted, bgzipped, tabix-indexed TSV.

    Args:
        csv_path: Path to input CSV file
        output_dir: Directory for output files (default: same as input)

    Returns:
        Path to the bgzipped file (.gz)

    Coordinate system: The output uses 1-based coordinates with columns:
        - Column 1: seqnames (chromosome)
        - Column 2: start (1-based)
        - Column 3: end (1-based, inclusive)
    """
    print(f"  Converting {csv_path} to tabix format...")

    # Read CSV
    df = pd.read_csv(csv_path)

    # Ensure numeric coordinates
    df['start'] = pd.to_numeric(df['start'], errors='coerce').astype('Int64')
    df['end'] = pd.to_numeric(df['end'], errors='coerce').astype('Int64')

    # Sort by chromosome and start position
    df['_chr_sort'] = df['seqnames'].apply(chr_sort_key)
    df = df.sort_values(['_chr_sort', 'start']).drop(columns=['_chr_sort'])

    # Determine output path
    if output_dir is None:
        output_dir = os.path.dirname(csv_path) or '.'

    base_name = os.path.splitext(os.path.basename(csv_path))[0]
    tsv_path = os.path.join(output_dir, f"{base_name}.sorted.tsv")
    gz_path = tsv_path + '.gz'

    # Ensure seqnames, start, end are first columns for tabix
    cols = ['seqnames', 'start', 'end'] + [c for c in df.columns if c not in ['seqnames', 'start', 'end']]
    df = df[cols]

    # Write TSV with header
    df.to_csv(tsv_path, sep='\t', index=False)

    # Bgzip compress
    print(f"  Compressing with bgzip...")
    try:
        subprocess.run(['bgzip', '-f', tsv_path], check=True, capture_output=True)
    except FileNotFoundError:
        raise RuntimeError(
            "bgzip not found. Please install htslib:\n"
            "  macOS: brew install htslib\n"
            "  Ubuntu: apt-get install tabix"
        )

    # Create tabix index
    print(f"  Creating tabix index...")
    try:
        # -s1: sequence column, -b2: begin column, -e3: end column, -S1: skip 1 header line
        subprocess.run(['tabix', '-s1', '-b2', '-e3', '-S1', '-f', gz_path], check=True, capture_output=True)
    except FileNotFoundError:
        raise RuntimeError(
            "tabix not found. Please install htslib:\n"
            "  macOS: brew install htslib\n"
            "  Ubuntu: apt-get install tabix"
        )

    print(f"  Created: {gz_path} and {gz_path}.tbi")
    return gz_path


def prepare_tabix_file(filepath: str, output_dir: Optional[str] = None) -> str:
    """
    Prepare a file for tabix queries. Converts CSV to tabix if needed.

    Args:
        filepath: Path to input file (CSV or .gz)
        output_dir: Directory for output files if conversion needed

    Returns:
        Path to tabix-ready .gz file
    """
    if is_tabix_indexed(filepath):
        print(f"  File already tabix-indexed: {filepath}")
        return filepath

    if filepath.endswith('.csv'):
        return csv_to_tabix(filepath, output_dir)

    if filepath.endswith('.gz'):
        # .gz but no .tbi - need to index it
        print(f"  Creating tabix index for {filepath}...")
        try:
            subprocess.run(['tabix', '-s1', '-b2', '-e3', '-S1', '-f', filepath], check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to index {filepath}. Ensure it's a sorted, bgzipped TSV.")
        return filepath

    raise ValueError(f"Unsupported file format: {filepath}. Expected .csv or .gz")


# =============================================================================
# Tabix Query Functions
# =============================================================================

def query_tabix_pysam(gz_path: str, chrom: str, start: int, end: int) -> list[dict]:
    """Query a tabix file using pysam."""
    results = []

    try:
        tbx = pysam.TabixFile(gz_path)

        # Try both with and without 'chr' prefix
        chroms_to_try = [chrom, normalize_chromosome(chrom), chrom.replace('chr', '')]
        chroms_to_try = list(set(chroms_to_try))  # Deduplicate

        header = None
        for line in tbx.header:
            header = line.split('\t')
            break

        for try_chrom in chroms_to_try:
            try:
                for row in tbx.fetch(try_chrom, start - 1, end):  # pysam uses 0-based
                    fields = row.split('\t')
                    if header:
                        results.append(dict(zip(header, fields)))
                    else:
                        results.append({'raw': fields})
                if results:
                    break
            except ValueError:
                continue  # Chromosome not in file

        tbx.close()
    except Exception as e:
        print(f"Warning: tabix query failed: {e}")

    return results


def query_tabix_subprocess(gz_path: str, chrom: str, start: int, end: int) -> list[dict]:
    """Query a tabix file using subprocess (fallback)."""
    results = []

    # Get header
    try:
        header_result = subprocess.run(
            ['tabix', '-H', gz_path],
            capture_output=True, text=True, check=True
        )
        header = header_result.stdout.strip().split('\t') if header_result.stdout.strip() else None
    except subprocess.CalledProcessError:
        header = None

    # Try both with and without 'chr' prefix
    chroms_to_try = [chrom, normalize_chromosome(chrom), chrom.replace('chr', '')]
    chroms_to_try = list(set(chroms_to_try))

    for try_chrom in chroms_to_try:
        region = f"{try_chrom}:{start}-{end}"
        try:
            result = subprocess.run(
                ['tabix', gz_path, region],
                capture_output=True, text=True, check=True
            )
            if result.stdout.strip():
                for line in result.stdout.strip().split('\n'):
                    fields = line.split('\t')
                    if header:
                        results.append(dict(zip(header, fields)))
                    else:
                        results.append({'raw': fields})
                if results:
                    break
        except subprocess.CalledProcessError:
            continue

    return results


def query_tabix(gz_path: str, chrom: str, start: int, end: int, padding: int = 1000) -> list[dict]:
    """
    Query a tabix-indexed file for records in a region.

    Args:
        gz_path: Path to bgzipped, tabix-indexed file
        chrom: Chromosome name
        start: Start position (1-based)
        end: End position (1-based, inclusive)
        padding: Padding to add around the region

    Returns:
        List of dictionaries with column names as keys
    """
    query_start = max(1, start - padding)
    query_end = end + padding

    if HAS_PYSAM:
        return query_tabix_pysam(gz_path, chrom, query_start, query_end)
    else:
        return query_tabix_subprocess(gz_path, chrom, query_start, query_end)


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_fraser_data(filepath: str) -> pd.DataFrame:
    """Load and prepare FRASER results data."""
    print(f"Loading FRASER data from {filepath}...")

    if filepath.endswith('.gz'):
        df = pd.read_csv(filepath, sep='\t', compression='gzip')
    else:
        df = pd.read_csv(filepath)

    validate_dataframe(df, FRASER_REQUIRED_COLUMNS, FRASER_OPTIONAL_COLUMNS, "FRASER")

    # Add hgncSymbol if missing
    if 'hgncSymbol' not in df.columns:
        df['hgncSymbol'] = 'NA'

    # Calculate -log10(pValue) for plotting
    df['-log10(pValue)'] = -np.log10(df['pValue'].replace(0, np.finfo(float).tiny))

    print(f"  Loaded {len(df)} FRASER rows")
    return df


def load_outrider_data(filepath: str) -> Optional[pd.DataFrame]:
    """Load and prepare OUTRIDER results data."""
    if filepath is None:
        return None

    print(f"Loading OUTRIDER data from {filepath}...")

    if filepath.endswith('.gz'):
        df = pd.read_csv(filepath, sep='\t', compression='gzip')
    else:
        df = pd.read_csv(filepath)

    validate_dataframe(df, OUTRIDER_REQUIRED_COLUMNS, OUTRIDER_OPTIONAL_COLUMNS, "OUTRIDER")

    # Standardize gene column name
    if 'hgncSymbol' not in df.columns and 'geneID' in df.columns:
        df['hgncSymbol'] = df['geneID']
    elif 'hgncSymbol' not in df.columns:
        df['hgncSymbol'] = 'NA'

    print(f"  Loaded {len(df)} OUTRIDER rows")
    return df


# =============================================================================
# Data Preparation for Dashboard
# =============================================================================

def get_top_positions(df: pd.DataFrame, n: int = 100) -> list[dict]:
    """Get top N positions by p-value for IGV."""
    top_df = df.nsmallest(n, 'pValue').copy()

    positions = []
    for _, row in top_df.iterrows():
        positions.append({
            'chr': str(row['seqnames']),
            'start': int(row['start']),
            'end': int(row['end']),
            'gene': str(row['hgncSymbol']) if pd.notna(row['hgncSymbol']) else 'NA',
            'pValue': float(row['pValue']),
            'deltaPsi': float(row['deltaPsi']),
            'psiValue': float(row['psiValue'])
        })

    return positions


def prepare_volcano_data(df: pd.DataFrame) -> dict:
    """Prepare data for volcano plot."""
    data = {
        'deltaPsi': df['deltaPsi'].tolist(),
        'log10pValue': df['-log10(pValue)'].tolist(),
        'padjust': df['padjust'].tolist(),
        'gene': df['hgncSymbol'].fillna('NA').tolist(),
        'chr': df['seqnames'].astype(str).tolist(),
        'start': df['start'].tolist(),
        'end': df['end'].tolist(),
        'psiValue': df['psiValue'].tolist(),
        'type': df['type'].tolist(),
        'sampleID': df['sampleID'].tolist(),
        'familyID': df['familyID'].fillna('Unknown').tolist() if 'familyID' in df.columns else ['Unknown'] * len(df)
    }

    # Add variant proximity annotations if available
    if 'closest_variant_pos' in df.columns:
        data['closest_variant_chr'] = df['closest_variant_chr'].fillna('').tolist()
        data['closest_variant_pos'] = df['closest_variant_pos'].fillna(0).astype(int).tolist()
        data['closest_variant_ref'] = df['closest_variant_ref'].fillna('').tolist()
        data['closest_variant_alt'] = df['closest_variant_alt'].fillna('').tolist()
        data['closest_variant_distance'] = df['closest_variant_distance'].fillna(-1).astype(int).tolist()

    return data


def prepare_outrider_volcano_data(df: pd.DataFrame) -> dict:
    """Prepare data for OUTRIDER volcano plot (zScore vs -log10(pValue))."""
    # Add -log10(pValue) column if not present
    if '-log10(pValue)' not in df.columns:
        df['-log10(pValue)'] = -np.log10(df['pValue'].replace(0, 1e-300))

    return {
        'zScore': df['zScore'].fillna(0).tolist() if 'zScore' in df.columns else [0] * len(df),
        'log10pValue': df['-log10(pValue)'].tolist(),
        'padjust': df['padjust'].tolist(),
        'gene': df['hgncSymbol'].fillna('NA').tolist() if 'hgncSymbol' in df.columns else df.get('geneID', pd.Series(['NA'] * len(df))).tolist(),
        'geneID': df['geneID'].tolist() if 'geneID' in df.columns else ['NA'] * len(df),
        'log2fc': df['l2fc'].fillna(0).tolist() if 'l2fc' in df.columns else (df['log2fc'].fillna(0).tolist() if 'log2fc' in df.columns else [0] * len(df)),
        'sampleID': df['sampleID'].tolist(),
        'familyID': df['familyID'].fillna('Unknown').tolist() if 'familyID' in df.columns else ['Unknown'] * len(df)
    }


def prepare_outrider_top_genes_data(df: pd.DataFrame, n: int = 50) -> dict:
    """
    Prepare data for OUTRIDER top genes plot.
    Shows the top N genes by significance.
    """
    # Get gene column (prefer hgncSymbol, fallback to geneID)
    if 'hgncSymbol' in df.columns:
        df = df.copy()
        df['gene_display'] = df['hgncSymbol'].fillna(df.get('geneID', 'NA'))
    else:
        df = df.copy()
        df['gene_display'] = df.get('geneID', 'NA')

    # Get top N genes by p-value
    top_df = df.nsmallest(n, 'pValue').copy()

    # Add -log10(pValue) if not present
    if '-log10(pValue)' not in top_df.columns:
        top_df['-log10(pValue)'] = -np.log10(top_df['pValue'].replace(0, 1e-300))

    return {
        'genes': top_df['gene_display'].tolist(),
        'log10pValue': top_df['-log10(pValue)'].tolist(),
        'padjust': top_df['padjust'].tolist(),
        'zScore': top_df['zScore'].fillna(0).tolist() if 'zScore' in top_df.columns else [0] * len(top_df),
        'log2fc': top_df['l2fc'].fillna(0).tolist() if 'l2fc' in top_df.columns else (top_df['log2fc'].fillna(0).tolist() if 'log2fc' in top_df.columns else [0] * len(top_df)),
        'sampleID': top_df['sampleID'].tolist(),
        'familyID': top_df['familyID'].fillna('Unknown').tolist() if 'familyID' in top_df.columns else ['Unknown'] * len(top_df)
    }


def prepare_manhattan_data(df: pd.DataFrame) -> tuple[dict, list, list]:
    """Prepare data for Manhattan plot with cumulative positions."""
    df = df.copy()

    # Sort by chromosome and position
    df['chr_num'] = df['seqnames'].apply(chr_sort_key)
    df_sorted = df.sort_values(['chr_num', 'start'])

    # Calculate cumulative positions
    df_sorted['middle'] = (df_sorted['start'] + df_sorted['end']) / 2
    chr_groups = df_sorted.groupby('seqnames', sort=False)

    cumulative_pos = []
    cumulative_length = 0
    chr_centers = {}

    for chr_name, group in chr_groups:
        chr_length = group['middle'].max() - group['middle'].min() + 1e6
        for idx in group.index:
            pos = df_sorted.loc[idx, 'middle'] - group['middle'].min() + cumulative_length
            cumulative_pos.append(pos)
        chr_centers[chr_name] = cumulative_length + chr_length / 2
        cumulative_length += chr_length

    manhattan_data = {
        'cumPos': cumulative_pos,
        'log10pValue': df_sorted['-log10(pValue)'].tolist(),
        'padjust': df_sorted['padjust'].tolist(),
        'gene': df_sorted['hgncSymbol'].fillna('NA').tolist(),
        'chr': df_sorted['seqnames'].astype(str).tolist(),
        'start': df_sorted['start'].tolist(),
        'end': df_sorted['end'].tolist(),
        'deltaPsi': df_sorted['deltaPsi'].tolist(),
        'psiValue': df_sorted['psiValue'].tolist(),
        'type': df_sorted['type'].tolist(),
        'sampleID': df_sorted['sampleID'].tolist(),
        'familyID': df_sorted['familyID'].fillna('Unknown').tolist() if 'familyID' in df_sorted.columns else ['Unknown'] * len(df_sorted)
    }

    # Add variant proximity annotations if available
    if 'closest_variant_pos' in df_sorted.columns:
        manhattan_data['closest_variant_chr'] = df_sorted['closest_variant_chr'].fillna('').tolist()
        manhattan_data['closest_variant_pos'] = df_sorted['closest_variant_pos'].fillna(0).astype(int).tolist()
        manhattan_data['closest_variant_ref'] = df_sorted['closest_variant_ref'].fillna('').tolist()
        manhattan_data['closest_variant_alt'] = df_sorted['closest_variant_alt'].fillna('').tolist()
        manhattan_data['closest_variant_distance'] = df_sorted['closest_variant_distance'].fillna(-1).astype(int).tolist()

    chr_labels = list(chr_centers.keys())
    chr_positions = list(chr_centers.values())

    return manhattan_data, chr_labels, chr_positions


def prepare_top_genes_data(df: pd.DataFrame, n: int = 100) -> list[dict]:
    """Prepare data for top genes plot."""
    top_df = df.nsmallest(n, 'pValue').copy()

    # Add -log10(pValue) if not present
    if '-log10(pValue)' not in top_df.columns:
        top_df['-log10(pValue)'] = -np.log10(top_df['pValue'].replace(0, 1e-300))

    return {
        'gene': top_df['hgncSymbol'].fillna('NA').tolist(),
        'log10pValue': top_df['-log10(pValue)'].tolist(),
        'padjust': top_df['padjust'].tolist(),
        'deltaPsi': top_df['deltaPsi'].tolist(),
        'psiValue': top_df['psiValue'].tolist(),
        'type': top_df['type'].tolist(),
        'sampleID': top_df['sampleID'].tolist(),
        'familyID': top_df['familyID'].fillna('Unknown').tolist() if 'familyID' in top_df.columns else ['Unknown'] * len(top_df)
    }


def prepare_tabix_data_for_js(df: pd.DataFrame, data_type: str = 'fraser') -> list[dict]:
    """
    Prepare data for JavaScript-side region queries.

    This creates a lightweight index that can be used for client-side
    region lookups without requiring server-side tabix queries.

    For OUTRIDER data, genomic coordinates (seqnames, start, end) are optional.
    Only records with coordinates will be included for region queries.
    """
    records = []

    for _, row in df.iterrows():
        # Check if genomic coordinates are present (required for region queries)
        if 'seqnames' not in row or pd.isna(row.get('seqnames')):
            continue
        if 'start' not in row or pd.isna(row.get('start')):
            continue
        if 'end' not in row or pd.isna(row.get('end')):
            continue

        record = {
            'chr': str(row['seqnames']),
            'start': int(row['start']),
            'end': int(row['end']),
            'gene': str(row.get('hgncSymbol', 'NA')) if pd.notna(row.get('hgncSymbol')) else 'NA',
            'pValue': float(row['pValue']) if pd.notna(row['pValue']) else None,
            'padjust': float(row['padjust']) if pd.notna(row['padjust']) else None,
        }

        if data_type == 'fraser':
            record['deltaPsi'] = float(row['deltaPsi']) if pd.notna(row.get('deltaPsi')) else None
            record['psiValue'] = float(row['psiValue']) if pd.notna(row.get('psiValue')) else None
            record['type'] = str(row['type']) if pd.notna(row.get('type')) else 'NA'

            # Add variant proximity information if available
            if 'closest_variant_pos' in row.index and pd.notna(row.get('closest_variant_pos')):
                record['closest_variant'] = {
                    'chr': str(row['closest_variant_chr']),
                    'pos': int(row['closest_variant_pos']),
                    'ref': str(row['closest_variant_ref']),
                    'alt': str(row['closest_variant_alt']),
                    'distance': int(row['closest_variant_distance'])
                }
        elif data_type == 'outrider':
            record['zScore'] = float(row['zScore']) if pd.notna(row.get('zScore')) else None
            record['log2fc'] = float(row['log2fc']) if pd.notna(row.get('log2fc')) else None

        records.append(record)

    return records


def calculate_stats(fraser_df: pd.DataFrame, outrider_df: Optional[pd.DataFrame],
                    pvalue_threshold: float) -> dict:
    """Calculate summary statistics for the dashboard."""
    stats = {
        'total_variants': len(fraser_df),
        'chromosomes': fraser_df['seqnames'].nunique(),
        'unique_genes': fraser_df['hgncSymbol'].nunique(),
        'significant_count': int((fraser_df['padjust'] < pvalue_threshold).sum())
    }

    if outrider_df is not None:
        stats['outrider_total'] = len(outrider_df)

    return stats


# =============================================================================
# HTML Generation
# =============================================================================

def render_dashboard(
    fraser_df: pd.DataFrame,
    outrider_df: Optional[pd.DataFrame],
    output_path: str,
    genome: str = 'hg38',
    pvalue_threshold: float = 0.05,
    deltapsi_threshold: float = 0.2,
    zscore_threshold: float = 2.0,
    bam_mapping: Optional[dict] = None
) -> None:
    """Render the dashboard HTML using Jinja2 template."""

    print("Preparing data for dashboard...")

    # Prepare all data structures
    volcano_data = prepare_volcano_data(fraser_df)
    manhattan_data, chr_labels, chr_positions = prepare_manhattan_data(fraser_df)
    top_positions = get_top_positions(fraser_df, n=100)

    # Prepare OUTRIDER plot data if available
    outrider_volcano_data = None
    outrider_top_genes_data = None
    if outrider_df is not None and len(outrider_df) > 0:
        outrider_volcano_data = prepare_outrider_volcano_data(outrider_df)
        outrider_top_genes_data = prepare_outrider_top_genes_data(outrider_df, n=50)

    # Prepare tabix-like data for JS queries
    fraser_tabix_data = prepare_tabix_data_for_js(fraser_df, 'fraser')
    outrider_tabix_data = prepare_tabix_data_for_js(outrider_df, 'outrider') if outrider_df is not None else []

    # Find common significant results if OUTRIDER data is available
    df_common = pd.DataFrame()
    if outrider_df is not None and len(outrider_df) > 0:
        df_common, _, _ = find_common_significant(fraser_df, outrider_df, pvalue_threshold, deltapsi_threshold)

    # Prepare table data
    print("Preparing table data...")

    # FRASER table
    fraser_table_cols = ['hgncSymbol', 'seqnames', 'start', 'end', 'type', 'pValue', 'padjust',
                         'psiValue', 'deltaPsi', 'sampleID']
    if 'familyID' in fraser_df.columns:
        fraser_table_cols.append('familyID')

    # Add variant proximity columns if available
    if 'closest_variant_pos' in fraser_df.columns:
        fraser_table_cols.extend(['closest_variant_chr', 'closest_variant_pos',
                                  'closest_variant_ref', 'closest_variant_alt',
                                  'closest_variant_distance'])

    fraser_table = prepare_table_data(fraser_df, fraser_table_cols)
    fraser_table_data = fraser_table.to_dict('records')
    fraser_table_columns = [{"data": col, "title": col} for col in fraser_table.columns]

    # OUTRIDER table
    outrider_table_data = []
    outrider_table_columns = []
    if outrider_df is not None and len(outrider_df) > 0:
        outrider_table_cols = ['sampleID', 'pValue', 'padjust']
        if 'geneID' in outrider_df.columns:
            outrider_table_cols.insert(0, 'geneID')
        elif 'hgncSymbol' in outrider_df.columns:
            outrider_table_cols.insert(0, 'hgncSymbol')
        if 'zScore' in outrider_df.columns:
            outrider_table_cols.append('zScore')
        if 'l2fc' in outrider_df.columns:
            outrider_table_cols.append('l2fc')
        elif 'log2fc' in outrider_df.columns:
            outrider_table_cols.append('log2fc')
        if 'familyID' in outrider_df.columns:
            outrider_table_cols.append('familyID')

        outrider_table = prepare_table_data(outrider_df, outrider_table_cols)
        outrider_table_data = outrider_table.to_dict('records')
        outrider_table_columns = [{"data": col, "title": col} for col in outrider_table.columns]

    # Common results table
    common_table_data = []
    common_table_columns = []
    if len(df_common) > 0:
        common_table = prepare_table_data(df_common)
        common_table_data = common_table.to_dict('records')
        common_table_columns = [{"data": col, "title": col} for col in common_table.columns]

    # Calculate stats
    stats = calculate_stats(fraser_df, outrider_df, pvalue_threshold)
    stats['common_count'] = len(df_common)

    # Get unique families for filtering (combine from both FRASER and OUTRIDER)
    families = set()
    if 'familyID' in fraser_df.columns:
        families.update(fraser_df['familyID'].unique().tolist())
    if outrider_df is not None and 'familyID' in outrider_df.columns:
        families.update(outrider_df['familyID'].unique().tolist())
    families = sorted([f for f in families if f != 'Unknown'])

    # Set up Jinja2 environment
    template_dir = Path(__file__).parent / 'templates'
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('interactive_dashboard.html.j2')

    # Render template
    print("Rendering HTML template...")
    html_content = template.render(
        volcano_data=volcano_data,
        manhattan_data=manhattan_data,
        outrider_volcano_data=outrider_volcano_data,
        outrider_top_genes_data=outrider_top_genes_data,
        top_positions=top_positions,
        chr_labels=chr_labels,
        chr_positions=chr_positions,
        fraser_tabix_data=fraser_tabix_data,
        outrider_tabix_data=outrider_tabix_data,
        fraser_table_data=fraser_table_data,
        fraser_table_columns=fraser_table_columns,
        outrider_table_data=outrider_table_data,
        outrider_table_columns=outrider_table_columns,
        common_table_data=common_table_data,
        common_table_columns=common_table_columns,
        stats=stats,
        families=families,
        genome=genome,
        default_pvalue_threshold=pvalue_threshold,
        default_deltapsi_threshold=deltapsi_threshold,
        default_zscore_threshold=zscore_threshold,
        sample_bam_mapping=bam_mapping or {}
    )

    # Write output
    with open(output_path, 'w') as f:
        f.write(html_content)

    print(f"  Dashboard written to: {output_path}")


# =============================================================================
# CLI Interface
# =============================================================================

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Create interactive FRASER/OUTRIDER results dashboard',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage with FRASER only
    python create_interactive_dashboard.py --fraser results.csv --output dashboard.html
    
    # With OUTRIDER results
    python create_interactive_dashboard.py --fraser fraser.csv --outrider outrider.csv --output dashboard.html
    
    # With variant proximity annotation
    python create_interactive_dashboard.py --fraser fraser.csv --variant-table variants.ht --output dashboard.html
    
    # Complete example with all features
    python create_interactive_dashboard.py --fraser fraser.csv --outrider outrider.csv \\
        --variant-table variants.ht --family-mapping families.csv --output dashboard.html
    
    # With pre-indexed tabix files
    python create_interactive_dashboard.py --fraser fraser.tsv.gz --outrider outrider.tsv.gz --output dashboard.html
    
    # Custom thresholds and genome
    python create_interactive_dashboard.py --fraser results.csv --output dashboard.html \\
        --genome hg19 --pvalue-threshold 0.01 --deltapsi-threshold 0.3
        """
    )

    parser.add_argument(
        '--fraser', '-f',
        required=True,
        help='Path to FRASER results file (CSV or tabix-indexed .gz)'
    )

    parser.add_argument(
        '--outrider', '-r',
        required=False,
        default=None,
        help='Path to OUTRIDER results file (CSV or tabix-indexed .gz)'
    )

    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output HTML file path'
    )

    parser.add_argument(
        '--genome', '-g',
        default='hg38',
        help='Reference genome for IGV (default: hg38)'
    )

    parser.add_argument(
        '--pvalue-threshold',
        type=float,
        default=0.05,
        help='Default p-value threshold (default: 0.05)'
    )

    parser.add_argument(
        '--deltapsi-threshold',
        type=float,
        default=0.2,
        help='Default |Delta PSI| threshold (default: 0.2)'
    )

    parser.add_argument(
        '--zscore-threshold',
        type=float,
        default=2.0,
        help='Default |Z-Score| threshold for OUTRIDER (default: 2.0)'
    )

    parser.add_argument(
        '--bam-mapping',
        default=None,
        help='Path to TSV file mapping sample IDs to BAM URLs (columns: sampleID, bam_url, bai_url)'
    )

    parser.add_argument(
        '--prepare-tabix',
        action='store_true',
        help='Convert CSV files to tabix-indexed format (keeps the indexed files)'
    )

    parser.add_argument(
        '--tabix-output-dir',
        default=None,
        help='Directory for tabix output files (default: same as input)'
    )

    parser.add_argument(
        '--family-mapping',
        default=None,
        help='Path to CPG-to-Family mapping CSV (rdnow-export-project-summary format)'
    )

    parser.add_argument(
        '--variant-table',
        default=None,
        help='Path to Hail table (.ht) with rare variants of interest for proximity annotation'
    )

    parser.add_argument(
        '--max-variant-distance',
        type=int,
        default=1000000,
        help='Maximum distance (bp) to consider for closest variant (default: 1000000)'
    )

    return parser.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_args()

    print("=" * 60)
    print("FRASER/OUTRIDER Interactive Dashboard Generator")
    print("=" * 60)

    # Load CPG to Family mapping if provided
    cpg_to_family = {}
    if args.family_mapping:
        cpg_to_family = load_cpg_to_family_mapping(args.family_mapping)

    # Load BAM mapping if provided
    bam_mapping = {}
    if args.bam_mapping:
        bam_mapping = load_bam_mapping(args.bam_mapping)

    # Optionally prepare tabix files
    fraser_path = args.fraser
    outrider_path = args.outrider

    if args.prepare_tabix:
        print("\nPreparing tabix-indexed files...")
        if args.fraser.endswith('.csv'):
            fraser_path = prepare_tabix_file(args.fraser, args.tabix_output_dir)
        if args.outrider and args.outrider.endswith('.csv'):
            outrider_path = prepare_tabix_file(args.outrider, args.tabix_output_dir)

    # Load data
    print("\nLoading data...")
    fraser_df = load_fraser_data(args.fraser)  # Use original for full data load
    outrider_df = load_outrider_data(args.outrider) if args.outrider else None

    # Add family IDs if mapping was provided
    if cpg_to_family:
        print("\nAdding family IDs...")
        fraser_df = add_family_ids(fraser_df, cpg_to_family)
        if outrider_df is not None:
            outrider_df = add_family_ids(outrider_df, cpg_to_family)

    # Load variant table and annotate with proximity if provided
    if args.variant_table:
        if not HAS_HAIL:
            print("\n⚠️  Warning: Hail is not installed. Skipping variant proximity annotation.")
            print("    Install Hail to enable this feature: pip install hail")
        else:
            try:
                variants_df = load_hail_variants(args.variant_table)
                fraser_df = find_closest_variant(
                    fraser_df,
                    variants_df,
                    max_distance=args.max_variant_distance
                )
            except Exception as e:
                print(f"\n⚠️  Warning: Failed to load variant table: {e}")
                print("    Continuing without variant proximity annotations.")

    # Get top positions info
    if len(fraser_df) > 0:
        top_positions = get_top_positions(fraser_df, n=1)
        if top_positions:
            print(f"\nTop FRASER hit: {top_positions[0]['chr']}:{top_positions[0]['start']}-{top_positions[0]['end']} "
                  f"(p={top_positions[0]['pValue']:.2e}, gene={top_positions[0]['gene']})")

    # Create dashboard
    print("\nCreating interactive dashboard...")
    render_dashboard(
        fraser_df=fraser_df,
        outrider_df=outrider_df,
        output_path=args.output,
        genome=args.genome,
        pvalue_threshold=args.pvalue_threshold,
        deltapsi_threshold=args.deltapsi_threshold,
        zscore_threshold=args.zscore_threshold,
        bam_mapping=bam_mapping
    )

    print("\n" + "=" * 60)
    print("✅ Dashboard created successfully!")
    print("=" * 60)
    print(f"\nOutput: {args.output}")
    print("\nFeatures:")
    print("  - Interactive threshold controls for p-value, |Delta PSI|, and |Z-Score|")
    print("  - Color-coded significance in volcano and Manhattan plots")
    print("  - Click points to view FRASER/OUTRIDER details and jump in IGV")
    print("  - IGV browser with top 100 most significant positions")
    print("  - Delta PSI track visualization in IGV")
    print("  - BAM loading controls for curators")
    print("  - Real-time statistics updates")
    if cpg_to_family:
        print(f"  - Family-based filtering enabled (searchable dropdown)")
    if bam_mapping:
        print(f"  - Pre-configured BAM files for {len(bam_mapping)} samples")
    if args.variant_table and 'closest_variant_pos' in fraser_df.columns:
        annotated_count = fraser_df['closest_variant_pos'].notna().sum()
        print(f"  - Variant proximity annotations ({annotated_count} Fraser results annotated)")

    if args.prepare_tabix:
        print("\nTabix files created:")
        if fraser_path != args.fraser:
            print(f"  - {fraser_path}")
        if outrider_path and outrider_path != args.outrider:
            print(f"  - {outrider_path}")


if __name__ == '__main__':
    main()

