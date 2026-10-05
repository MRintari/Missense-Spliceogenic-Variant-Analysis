# Import Modules
import pandas as pd
import numpy as np
import pysam

import pyranges as pr

# Load ClinVar Files
# Created on command line bash filtering for each CLNSIG
benign = pysam.VariantFile('/slade/home/mr935/data/clinvar_annots/benign_clinvar_all_vars_indelslt20bp.bcf')
likely_benign = pysam.VariantFile('/slade/home/mr935/data/clinvar_annots/likely_benign_clinvar_all_vars_indelslt20bp.bcf')
vus = pysam.VariantFile('/slade/home/mr935/data/clinvar_annots/vus_clinvar_all_vars_indelslt20bp.bcf')
likely_pathogenic = pysam.VariantFile('/slade/home/mr935/data/clinvar_annots/likely_pathogenic_clinvar_all_vars_indelslt20bp.bcf')
pathogenic = pysam.VariantFile('/slade/home/mr935/data/clinvar_annots/pathogenic_clinvar_all_vars_indelslt20bp.bcf')

# Identify correct transcript for each gene symbol to have one set of SpliceAI annotations per ID
mane = pd.read_csv('/slade/home/mr935/data/MANE.GRCh38.v1.3.summary.txt.gz', sep='\t')
mane_data = mane[mane['MANE_status'] == 'MANE Select'].copy()
mane_data['MANE'] = mane_data['Ensembl_nuc'].str.split('.').str[0]

def extract_missense_gene(info):
    # Get the CSQ tuple from the dictionary
    csq_data = info.get('CSQ', [])
    
    # Find the first string containing 'missense'
    match = next((s for s in csq_data if 'missense' in s), None)
    
    # If a match exists, split and return the 4th element (index 3)
    if match:
        return match.split('|')[3]
    
    return None # Or "Unknown" if no missense is found


# Clean & Filter ClinVar 
def clean_clinvar(df, mane_data=mane_data):
    """
    Pipeline to clean ClinVar data and extract their
    SpliceAI scores. Handles multi-trancript genes, 
    conflicting predicted coding consequence, missing scores,
    etc.
    """
    # Create copy to avoid modifying original data
    df = df.copy()
    # Ensure data is in string format
    df['spliceai'] = df['spliceai'].astype(str)

    # Remove leading and trailing brackets
    df['spliceai'] = df['spliceai'].apply(
        lambda x: x[1:-1] if isinstance(x, (tuple, list, str)) else None
    )
    
    # Extract the appropriate gene for the variant
    df['gene'] = df['info'].apply(extract_missense_gene)

    # Extract the correct transcript for the variant's gene
    df = df.merge(
        mane_data[['symbol', 'MANE']],
        left_on='gene',
        right_on='symbol',
        how='left'
    )

    # Filter for missense variants
    df_missense = df[df['info'].apply(lambda x: 'missense_variant' in str(x.values()))].copy()

    # Filter variants in which MC specifies missense but CSQ does not state it as missense
    df_missense = df_missense[df_missense['gene'].notna()].copy()
    df_missense['spliceai'] = df_missense['spliceai'].astype('str')

    # Identify variants in multiple genes
    multi_tx = df_missense[df_missense['spliceai'].str.count(r'\|') > 20].copy()

    # Identify correct set of SpliceAI scores for the variant
    multi_tx['spliceai'] = multi_tx.apply(extract_matching_spliceai, axis=1)

    # Add cleaned data to original df
    df_missense = df_missense[~(df_missense['id'].isin(multi_tx['id']))].copy()
    clean_df_missense = pd.concat([df_missense, multi_tx], ignore_index=True)

    # Expand SpliceAI vals to columns
    splice_cols = ['ALT', 'Transcript', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL', 'RS_AG',
               'RS_AL', 'RS_DG', 'RS_DL', 'MANEselectDonorSpliceSitesWithinContext.DP.RS_REF.RS_ALT.',
               'MANEselectAcceptorSpliceSitesWithincontext.DP.RS_REF.RS_ALT.', 'DonorSitesWithRawScoregt0.5.DP.RS_REF.RS_ALT.',
               'AcceptorSitesWithRawScoregt0.5.DP.RS_REF.RS_ALT.', 'SPLICE_CHANGE.Donor.Acceptor.',
               'FRAME_CHANGE.Donor.Acceptor.',	'AA_CHANGE.Donor.Acceptor.']
    
    splice_info = clean_df_missense['spliceai'].str.split('|', expand=True)
    splice_info.columns = splice_cols
    clean_df_missense = pd.concat([clean_df_missense.drop(columns=['spliceai']), splice_info], axis=1)

    return clean_df_missense

clean_benign_missense = clean_clinvar(benign)
clean_lb_missense = clean_clinvar(likely_benign)
clean_vus_missense = clean_clinvar(vus)
clean_lp_missense = clean_clinvar(likely_pathogenic)
clean_pathogenic_missense = clean_clinvar(pathogenic)

# Save for later use
clean_benign_missense.to_parquet('/slade/home/mr935/data/clinvar_annots/clean_missense_annots/clean_benign_missense.parquet')
clean_lb_missense.to_parquet('/slade/home/mr935/data/clinvar_annots/clean_missense_annots/clean_lb_missense.parquet')
clean_vus_missense.to_parquet('/slade/home/mr935/data/clinvar_annots/clean_missense_annots/clean_vus_missense.parquet')
clean_lp_missense.to_parquet('/slade/home/mr935/data/clinvar_annots/clean_missense_annots/clean_lp_missense.parquet')
clean_pathogenic_missense.to_parquet('/slade/home/mr935/data/clinvar_annots/clean_missense_annots/clean_pathogenic_missense.parquet')
