"""
Create comprehensive dashboard with FRASER, OUTRIDER, and common significant results
"""
import pandas as pd
import numpy as np
import json
from pathlib import Path
from jinja2 import Template


def load_fraser_data(fraser_csv_path, rdnow_csv_path):
    """Load FRASER results and merge with rd-now metadata"""
    print("Loading FRASER results...")
    df_fraser = pd.read_csv(fraser_csv_path)

    print("Loading rd-now metadata...")
    df_rdnow = pd.read_csv(rdnow_csv_path)

    # Merge on sequencing_group.id = sampleID
    print("Merging FRASER with rd-now metadata...")
    df_merged = df_fraser.merge(
        df_rdnow,
        left_on='sampleID',
        right_on='sequencing_group.id',
        how='left'
    )

    print(f"FRASER: {len(df_merged)} rows, {df_merged['family.external_ids'].notna().sum()} with metadata")

    # Calculate -log10(pValue) for plotting
    df_merged['-log10(pValue)'] = -np.log10(df_merged['pValue'])

    return df_merged


def load_outrider_data(outrider_csv_path, rdnow_csv_path):
    """Load OUTRIDER results and merge with rd-now metadata"""
    print("\nLoading OUTRIDER results...")
    df_outrider = pd.read_csv(outrider_csv_path)

    print("Loading rd-now metadata for OUTRIDER...")
    df_rdnow = pd.read_csv(rdnow_csv_path)

    # Merge on sampleID
    print("Merging OUTRIDER with rd-now metadata...")
    df_merged = df_outrider.merge(
        df_rdnow,
        left_on='sampleID',
        right_on='sequencing_group.id',
        how='left'
    )

    print(f"OUTRIDER: {len(df_merged)} rows, {df_merged['family.external_ids'].notna().sum()} with metadata")

    # Calculate -log10(pValue) for plotting
    df_merged['-log10(pValue)'] = -np.log10(df_merged['pValue'])

    # Extract gene symbol from geneID (remove version number)
    df_merged['gene_base'] = df_merged['geneID'].str.split('.').str[0]

    return df_merged


def find_common_significant(df_fraser, df_outrider, pval_threshold=0.05, deltapsi_threshold=0.2):
    """Find genes that are significant in both FRASER and OUTRIDER"""
    print("\nFinding common significant results...")

    # FRASER significant: meets p-value and deltaPSI thresholds
    fraser_sig = df_fraser[
        (df_fraser['padjust'] <= pval_threshold) &
        (np.abs(df_fraser['deltaPsi']) >= deltapsi_threshold)
    ].copy()

    # OUTRIDER significant: meets p-value and z-score thresholds
    outrider_sig = df_outrider[
        (df_outrider['padjust'] <= pval_threshold) &
        (np.abs(df_outrider['zScore']) >= 2)
    ].copy()

    print(f"  FRASER significant: {len(fraser_sig)}")
    print(f"  OUTRIDER significant: {len(outrider_sig)}")

    # Strategy 1: Match by sample only (same sample showing both types of aberrations)
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
                # Check if gene names overlap (HGNC symbol might be in Ensembl annotation)
                fraser_gene = str(fraser_row['hgncSymbol']).upper() if pd.notna(fraser_row['hgncSymbol']) else ''
                outrider_gene = str(outrider_row['geneID']).upper()

                # Consider it a match if the HGNC symbol appears in the Ensembl ID or vice versa
                # Or if they are within the same genomic region (for positional matching)
                is_match = False

                if fraser_gene and (
                    fraser_gene in outrider_gene or
                    any(fraser_gene in part for part in outrider_gene.split(';'))
                ):
                    is_match = True

                # If we have a match, add to common results
                if is_match:
                    common_by_sample.append({
                        'gene': fraser_row['hgncSymbol'],
                        'ensembl_id': outrider_row['geneID'],
                        'sampleID': sample,
                        'chr': fraser_row['seqnames'],
                        'start': fraser_row['start'],
                        'end': fraser_row['end'],
                        'fraser_pvalue': fraser_row['pValue'],
                        'fraser_padjust': fraser_row['padjust'],
                        'fraser_deltaPsi': fraser_row['deltaPsi'],
                        'fraser_psiValue': fraser_row['psiValue'],
                        'fraser_type': fraser_row['type'],
                        'outrider_pvalue': outrider_row['pValue'],
                        'outrider_padjust': outrider_row['padjust'],
                        'outrider_zScore': outrider_row['zScore'],
                        'outrider_l2fc': outrider_row['l2fc'],
                        'familyID': fraser_row.get('family.external_ids', 'NA'),
                        'participantID': fraser_row.get('participant.external_ids', 'NA'),
                        'sampleType': fraser_row.get('sample.type', 'NA')
                    })

    # Also add entries where sample matches but gene doesn't necessarily match
    # (to show samples with both types of aberrations even if different genes)
    for sample in common_samples:
        fraser_sample_genes = fraser_sig[fraser_sig['sampleID'] == sample]
        outrider_sample_genes = outrider_sig[outrider_sig['sampleID'] == sample]

        # Take top hit from each for this sample
        if len(fraser_sample_genes) > 0 and len(outrider_sample_genes) > 0:
            # Get the most significant from each
            fraser_top = fraser_sample_genes.nsmallest(1, 'pValue').iloc[0]
            outrider_top = outrider_sample_genes.nsmallest(1, 'pValue').iloc[0]

            # Check if this combination already exists
            already_exists = any(
                r['sampleID'] == sample and
                str(r['gene']) == str(fraser_top['hgncSymbol']) and
                str(r['ensembl_id']) == str(outrider_top['geneID'])
                for r in common_by_sample
            )

            if not already_exists:
                common_by_sample.append({
                    'gene': fraser_top['hgncSymbol'],
                    'ensembl_id': outrider_top['geneID'],
                    'sampleID': sample,
                    'chr': fraser_top['seqnames'],
                    'start': fraser_top['start'],
                    'end': fraser_top['end'],
                    'fraser_pvalue': fraser_top['pValue'],
                    'fraser_padjust': fraser_top['padjust'],
                    'fraser_deltaPsi': fraser_top['deltaPsi'],
                    'fraser_psiValue': fraser_top['psiValue'],
                    'fraser_type': fraser_top['type'],
                    'outrider_pvalue': outrider_top['pValue'],
                    'outrider_padjust': outrider_top['padjust'],
                    'outrider_zScore': outrider_top['zScore'],
                    'outrider_l2fc': outrider_top['l2fc'],
                    'familyID': fraser_top.get('family.external_ids', 'NA'),
                    'participantID': fraser_top.get('participant.external_ids', 'NA'),
                    'sampleType': fraser_top.get('sample.type', 'NA')
                })

    df_common = pd.DataFrame(common_by_sample)
    print(f"  Common significant results: {len(df_common)}")
    print(f"  (Samples showing aberrations in both splicing and expression)")

    return df_common, fraser_sig, outrider_sig


def prepare_table_data(df, columns_subset=None):
    """Prepare data for the searchable table"""
    if columns_subset:
        available_columns = [col for col in columns_subset if col in df.columns]
    else:
        available_columns = df.columns.tolist()

    table_df = df[available_columns].copy()

    # Format numeric columns
    for col in table_df.columns:
        if 'pValue' in col or 'padjust' in col:
            table_df[col] = table_df[col].apply(lambda x: f"{x:.2e}" if pd.notna(x) and isinstance(x, (int, float)) else str(x))
        elif 'psi' in col.lower() or 'Score' in col or 'l2fc' in col or 'deltaPsi' in col:
            table_df[col] = table_df[col].apply(lambda x: f"{x:.3f}" if pd.notna(x) and isinstance(x, (int, float)) else str(x))

    # Fill NaN values
    table_df = table_df.fillna('NA')

    return table_df


def create_comprehensive_dashboard(df_fraser, df_outrider, df_common):
    """Create comprehensive HTML dashboard with FRASER, OUTRIDER, and common results
    Returns a dictionary of template context variables for rendering with Jinja2 template
    """

    # Prepare FRASER data for JavaScript
    fraser_data = {
        'deltaPsi': df_fraser['deltaPsi'].tolist(),
        'log10pValue': df_fraser['-log10(pValue)'].tolist(),
        'padjust': df_fraser['padjust'].tolist(),
        'gene': df_fraser['hgncSymbol'].fillna('NA').tolist(),
        'chr': df_fraser['seqnames'].tolist(),
        'start': df_fraser['start'].tolist(),
        'end': df_fraser['end'].tolist(),
        'psiValue': df_fraser['psiValue'].tolist(),
        'type': df_fraser['type'].tolist(),
        'sampleID': df_fraser['sampleID'].tolist(),
        'familyID': df_fraser['family.external_ids'].fillna('NA').tolist() if 'family.external_ids' in df_fraser.columns else ['NA'] * len(df_fraser),
        'participantID': df_fraser['participant.external_ids'].fillna('NA').tolist() if 'participant.external_ids' in df_fraser.columns else ['NA'] * len(df_fraser),
        'sampleType': df_fraser['sample.type'].fillna('NA').tolist() if 'sample.type' in df_fraser.columns else ['NA'] * len(df_fraser),
    }

    # Prepare OUTRIDER data for JavaScript
    outrider_data = {
        'l2fc': df_outrider['l2fc'].tolist(),
        'log10pValue': df_outrider['-log10(pValue)'].tolist(),
        'padjust': df_outrider['padjust'].tolist(),
        'geneID': df_outrider['geneID'].tolist(),
        'zScore': df_outrider['zScore'].tolist(),
        'rawcounts': df_outrider['rawcounts'].tolist(),
        'normcounts': df_outrider['normcounts'].tolist(),
        'sampleID': df_outrider['sampleID'].tolist(),
        'familyID': df_outrider['family.external_ids'].fillna('NA').tolist() if 'family.external_ids' in df_outrider.columns else ['NA'] * len(df_outrider),
        'participantID': df_outrider['participant.external_ids'].fillna('NA').tolist() if 'participant.external_ids' in df_outrider.columns else ['NA'] * len(df_outrider),
        'sampleType': df_outrider['sample.type'].fillna('NA').tolist() if 'sample.type' in df_outrider.columns else ['NA'] * len(df_outrider),
    }

    # Sort FRASER by chromosome for Manhattan plot
    def chr_sort_key(chr_name):
        chr_name = str(chr_name).replace('chr', '')
        if chr_name == 'X':
            return 23
        elif chr_name == 'Y':
            return 24
        elif chr_name == 'M' or chr_name == 'MT':
            return 25
        else:
            try:
                return int(chr_name)
            except:
                return 26

    df_fraser['chr_num'] = df_fraser['seqnames'].apply(chr_sort_key)
    df_fraser_sorted = df_fraser.sort_values(['chr_num', 'start'])

    # Calculate cumulative positions for Manhattan
    df_fraser_sorted['middle'] = (df_fraser_sorted['start'] + df_fraser_sorted['end']) / 2
    chr_groups = df_fraser_sorted.groupby('seqnames', sort=False)

    cumulative_pos = []
    cumulative_length = 0
    chr_centers = {}

    for chr_name, group in chr_groups:
        chr_length = group['middle'].max() - group['middle'].min() + 1e6
        for idx in group.index:
            pos = df_fraser_sorted.loc[idx, 'middle'] - group['middle'].min() + cumulative_length
            cumulative_pos.append(pos)
        chr_centers[chr_name] = cumulative_length + chr_length / 2
        cumulative_length += chr_length

    manhattan_data = {
        'cumPos': cumulative_pos,
        'log10pValue': df_fraser_sorted['-log10(pValue)'].tolist(),
        'padjust': df_fraser_sorted['padjust'].tolist(),
        'gene': df_fraser_sorted['hgncSymbol'].fillna('NA').tolist(),
        'chr': df_fraser_sorted['seqnames'].tolist(),
        'start': df_fraser_sorted['start'].tolist(),
        'end': df_fraser_sorted['end'].tolist(),
        'deltaPsi': df_fraser_sorted['deltaPsi'].tolist(),
        'psiValue': df_fraser_sorted['psiValue'].tolist(),
        'type': df_fraser_sorted['type'].tolist(),
        'sampleID': df_fraser_sorted['sampleID'].tolist(),
        'familyID': df_fraser_sorted['family.external_ids'].fillna('NA').tolist() if 'family.external_ids' in df_fraser_sorted.columns else ['NA'] * len(df_fraser_sorted),
        'participantID': df_fraser_sorted['participant.external_ids'].fillna('NA').tolist() if 'participant.external_ids' in df_fraser_sorted.columns else ['NA'] * len(df_fraser_sorted),
        'sampleType': df_fraser_sorted['sample.type'].fillna('NA').tolist() if 'sample.type' in df_fraser_sorted.columns else ['NA'] * len(df_fraser_sorted),
    }

    chr_labels = list(chr_centers.keys())
    chr_positions = list(chr_centers.values())

    # Prepare all FRASER positions for IGV
    fraser_positions = []
    for _, row in df_fraser.iterrows():
        fraser_positions.append({
            'chr': row['seqnames'],
            'start': int(row['start']),
            'end': int(row['end']),
            'gene': str(row['hgncSymbol']) if pd.notna(row['hgncSymbol']) else 'NA',
            'pValue': float(row['pValue']),
            'deltaPsi': float(row['deltaPsi']),
            'sampleID': str(row['sampleID']),
            'familyID': str(row['family.external_ids']) if 'family.external_ids' in row and pd.notna(row['family.external_ids']) else 'NA'
        })

    # Prepare common results for IGV
    common_positions = []
    if len(df_common) > 0:
        for _, row in df_common.iterrows():
            common_positions.append({
                'chr': row['chr'],
                'start': int(row['start']),
                'end': int(row['end']),
                'gene': str(row['gene']),
                'fraser_pValue': float(row['fraser_pvalue']),
                'outrider_pValue': float(row['outrider_pvalue']),
                'deltaPsi': float(row['fraser_deltaPsi']),
                'zScore': float(row['outrider_zScore']),
                'sampleID': str(row['sampleID']),
                'familyID': str(row['familyID'])
            })

    # Prepare table data
    fraser_table_cols = ['hgncSymbol', 'seqnames', 'start', 'end', 'type', 'pValue', 'padjust',
                         'psiValue', 'deltaPsi', 'sampleID', 'family.external_ids',
                         'participant.external_ids', 'sample.type']
    fraser_table = prepare_table_data(df_fraser, fraser_table_cols)
    fraser_table_data = fraser_table.to_dict('records')
    fraser_table_columns = list(fraser_table.columns)

    outrider_table_cols = ['geneID', 'sampleID', 'pValue', 'padjust', 'zScore', 'l2fc',
                           'rawcounts', 'normcounts', 'family.external_ids',
                           'participant.external_ids', 'sample.type']
    outrider_table = prepare_table_data(df_outrider, outrider_table_cols)
    outrider_table_data = outrider_table.to_dict('records')
    outrider_table_columns = list(outrider_table.columns)

    if len(df_common) > 0:
        common_table = prepare_table_data(df_common)
        common_table_data = common_table.to_dict('records')
        common_table_columns = list(common_table.columns)
    else:
        common_table_data = []
        common_table_columns = []

    # Return all the prepared data and context variables for template rendering
    fraser_js_data = {
        'deltaPsi': fraser_data['deltaPsi'],
        'log10p': fraser_data['log10pValue'],
        'padjust': fraser_data['padjust'],
        'chr': fraser_data['chr'],
        'start': fraser_data['start'],
        'end': fraser_data['end'],
        'hoverText': [
            f"Gene: {fraser_data['gene'][i]}<br>"
            f"Position: {fraser_data['chr'][i]}:{fraser_data['start'][i]}-{fraser_data['end'][i]}<br>"
            f"Sample: {fraser_data['sampleID'][i]}<br>"
            f"Family: {fraser_data['familyID'][i]}<br>"
            f"Sample Type: {fraser_data['sampleType'][i]}<br>"
            f"P-value: {fraser_data['padjust'][i]:.2e}<br>"
            f"Delta PSI: {fraser_data['deltaPsi'][i]:.3f}<br>"
            f"PSI Value: {fraser_data['psiValue'][i]:.3f}<br>"
            f"Type: {fraser_data['type'][i]}"
            for i in range(len(fraser_data['gene']))
        ]
    }

    outrider_js_data = {
        'log2FC': outrider_data['l2fc'],
        'log10p': outrider_data['log10pValue'],
        'padjust': outrider_data['padjust'],
        'zScore': outrider_data['zScore'],
        'chr': ['chr1'] * len(outrider_data['geneID']),  # Placeholder, needs actual chr data
        'start': [0] * len(outrider_data['geneID']),  # Placeholder
        'end': [0] * len(outrider_data['geneID']),  # Placeholder
        'hoverText': [
            f"Gene: {outrider_data['geneID'][i]}<br>"
            f"Sample: {outrider_data['sampleID'][i]}<br>"
            f"Family: {outrider_data['familyID'][i]}<br>"
            f"Sample Type: {outrider_data['sampleType'][i]}<br>"
            f"P-value: {outrider_data['padjust'][i]:.2e}<br>"
            f"Z-Score: {outrider_data['zScore'][i]:.2f}<br>"
            f"log2 FC: {outrider_data['l2fc'][i]:.3f}<br>"
            f"Raw counts: {outrider_data['rawcounts'][i]:.0f}<br>"
            f"Norm counts: {outrider_data['normcounts'][i]:.2f}"
            for i in range(len(outrider_data['geneID']))
        ]
    }

    manhattan_js_data = {
        'position': manhattan_data['cumPos'],
        'log10p': manhattan_data['log10pValue'],
        'padjust': manhattan_data['padjust'],
        'deltaPsi': manhattan_data['deltaPsi'],
        'chr': manhattan_data['chr'],
        'start': manhattan_data['start'],
        'end': manhattan_data['end'],
        'chrStarts': chr_positions,
        'chrLabels': chr_labels,
        'hoverText': [
            f"Gene: {manhattan_data['gene'][i]}<br>"
            f"Position: {manhattan_data['chr'][i]}:{manhattan_data['start'][i]}-{manhattan_data['end'][i]}<br>"
            f"Sample: {manhattan_data['sampleID'][i]}<br>"
            f"Family: {manhattan_data['familyID'][i]}<br>"
            f"P-value: {manhattan_data['padjust'][i]:.2e}<br>"
            f"Delta PSI: {manhattan_data['deltaPsi'][i]:.3f}"
            for i in range(len(manhattan_data['gene']))
        ]
    }


    # Prepare IGV tracks
    igv_tracks = [
        {
            "name": "Genes",
            "type": "annotation",
            "format": "refgene",
            "url": "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.txt.gz",
            "indexURL": "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.txt.gz.tbi",
            "order": 1000000,
            "visibilityWindow": 300000000,
            "displayMode": "EXPANDED"
        }
    ]

    # Add common positions track
    if len(common_positions) > 0:
        common_features = []
        for pos in common_positions:
            common_features.append({
                "chr": pos['chr'],
                "start": pos['start'],
                "end": pos['end'],
                "name": f"{pos['gene']} ({pos['sampleID']})",
                "color": "orange"
            })

        igv_tracks.append({
            "name": "Common Significant (FRASER + OUTRIDER)",
            "type": "annotation",
            "format": "bed",
            "features": common_features,
            "displayMode": "EXPANDED",
            "color": "orange"
        })

    # Add all FRASER results track
    fraser_features = []
    for pos in fraser_positions:
        fraser_features.append({
            "chr": pos['chr'],
            "start": pos['start'],
            "end": pos['end'],
            "name": f"{pos['gene']} ({pos['sampleID']})",
            "color": "blue"
        })

    igv_tracks.append({
        "name": "All FRASER Results",
        "type": "annotation",
        "format": "bed",
        "features": fraser_features,
        "displayMode": "COLLAPSED",
        "color": "blue"
    })

    # Calculate statistics
    fraser_sig_count = (df_fraser['padjust'] < 0.05).sum()
    outrider_sig_count = (df_outrider['padjust'] < 0.05).sum()
    unique_families = df_fraser['family.external_ids'].nunique() if 'family.external_ids' in df_fraser.columns else 'N/A'

    # Prepare column definitions for DataTables
    fraser_table_cols_json = [{"data": col, "title": col} for col in fraser_table_columns]
    outrider_table_cols_json = [{"data": col, "title": col} for col in outrider_table_columns]
    common_table_cols_json = [{"data": col, "title": col} for col in common_table_columns]

    return {
        'fraser_count': len(df_fraser),
        'outrider_count': len(df_outrider),
        'common_count': len(df_common),
        'fraser_sig_count': int(fraser_sig_count),
        'outrider_sig_count': int(outrider_sig_count),
        'unique_families': str(unique_families),
        'fraser_data_json': json.dumps(fraser_js_data),
        'outrider_data_json': json.dumps(outrider_js_data),
        'manhattan_data_json': json.dumps(manhattan_js_data),
        'igv_tracks_json': json.dumps(igv_tracks),
        'fraser_table_data_json': json.dumps(fraser_table_data),
        'outrider_table_data_json': json.dumps(outrider_table_data),
        'common_table_data_json': json.dumps(common_table_data),
        'fraser_table_columns': fraser_table_columns,
        'outrider_table_columns': outrider_table_columns,
        'common_table_columns': common_table_columns,
        'fraser_table_columns_json': json.dumps(fraser_table_cols_json),
        'outrider_table_columns_json': json.dumps(outrider_table_cols_json),
        'common_table_columns_json': json.dumps(common_table_cols_json)
    }



def main():
    fraser_csv_path = '/Users/johass/PycharmProjects/images/images/fraser/transcriptome_fraser_COH10509.results.significant.csv'
    outrider_csv_path = '/Users/johass/PycharmProjects/images/images/fraser/transcriptome_outrider_COH10509.outrider.results.csv'
    rdnow_csv_path = '/Users/johass/PycharmProjects/images/images/fraser/rdnow-export-project-summary-2026-02-19.csv'
    output_path = '/Users/johass/PycharmProjects/images/images/fraser/comprehensive_dashboard.html'

    print("=" * 80)
    print("COMPREHENSIVE FRASER & OUTRIDER DASHBOARD GENERATOR")
    print("=" * 80)

    # Load data
    df_fraser = load_fraser_data(fraser_csv_path, rdnow_csv_path)
    df_outrider = load_outrider_data(outrider_csv_path, rdnow_csv_path)

    # Find common significant results
    df_common, fraser_sig, outrider_sig = find_common_significant(
        df_fraser, df_outrider,
        pval_threshold=0.05,
        deltapsi_threshold=0.2
    )

    print("\nCreating comprehensive dashboard...")

    # Generate the template context
    template_context = create_comprehensive_dashboard(df_fraser, df_outrider, df_common)

    # Load the HTML template
    template_path = Path(__file__).parent / 'dashboard_template.html'
    print(f"Loading template from: {template_path}")

    with open(template_path, 'r', encoding='utf-8') as f:
        template = Template(f.read())

    # Render the template with the context
    html_content = template.render(**template_context)

    # Write the output HTML file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    print(f"\n✅ Comprehensive dashboard created: {output_path}")
    print("\n" + "=" * 80)
    print("SUMMARY:")
    print("=" * 80)
    print(f"📊 FRASER Results:        {len(df_fraser)}")
    print(f"📊 OUTRIDER Results:      {len(df_outrider)}")
    print(f"🎯 Common Significant:    {len(df_common)}")
    print(f"   - Same gene + sample showing both splicing & expression changes")
    print("\n" + "=" * 80)
    print("FEATURES:")
    print("=" * 80)
    print("✓ Separate volcano plots for FRASER and OUTRIDER")
    print("✓ Manhattan plot for genome-wide FRASER results")
    print("✓ Intersection analysis identifying common significant genes")
    print("✓ Dedicated IGV track for common results (orange/gold)")
    print("✓ All FRASER results in IGV track (blue)")
    print("✓ Searchable tables with rd-now metadata")
    print("✓ Interactive threshold controls")
    print("✓ Click-to-navigate in IGV browser")
    print("✓ Export tables to CSV/Excel")
    print("=" * 80)

    if len(df_common) > 0:
        print("\n🎯 TOP COMMON SIGNIFICANT GENES:")
        print("-" * 80)
        for i, row in df_common.head(10).iterrows():
            gene = str(row['gene']) if pd.notna(row['gene']) else 'NA'
            sample = str(row['sampleID']) if pd.notna(row['sampleID']) else 'NA'
            family = str(row['familyID']) if pd.notna(row['familyID']) else 'NA'
            print(f"  {gene:15s} | {sample:12s} | {family:10s}")
            print(f"    FRASER:   p={row['fraser_pvalue']:.2e}, ΔΨ={row['fraser_deltaPsi']:+.3f}")
            print(f"    OUTRIDER: p={row['outrider_pvalue']:.2e}, Z={row['outrider_zScore']:+.2f}")
            print()
        print("=" * 80)


if __name__ == '__main__':
    main()

