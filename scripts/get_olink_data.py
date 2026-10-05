# Import Modules
import pandas as pd

####################################################
####### CONFIG & HELPER FUNCTIONS ##################
####################################################

DATA_DIR = '/slade/home/mr935/data'
BIOBANK_DIR = f'{DATA_DIR}/biobank_data'
OLINK_DIR = f'{BIOBANK_DIR}/olink_files'


def format_and_export_mask(df, gene_col, var_col, output_path):
    """
    Formats gene symbols to lowercase, normalizes variant IDs (e.g. 1:1234:A:G),
    sorts by gene, and saves to TSV for extract_olink_data.sh.
    """
    mask = df[[gene_col, var_col]].copy()
    mask.columns = ['SYMBOL', 'ID']
    mask['SYMBOL'] = mask['SYMBOL'].str.lower()
    mask['ID'] = (
        mask['ID']
        .astype(str)
        .str.replace('chr', '', regex=False)
        .str.replace('-', ':', regex=False)
    )
    mask = mask.sort_values('SYMBOL').reset_index(drop=True)
    mask.to_csv(output_path, sep='\t', index=False)
    return mask


def clean_olink_output(file_path, mask_label):
    """
    Loads TSV output from extract_olink_data.sh, filters non-null Betas,
    standardizes gene/symbol column naming and upper-case genes, and assigns mask category.
    """
    df = pd.read_csv(file_path, sep='\t', index_col=False)
    df = df[~df['Beta'].isna()].reset_index(drop=True)

    # Standardize column naming for Gene
    if 'SYMBOL' in df.columns:
        df.rename(columns={'SYMBOL': 'Gene'}, inplace=True)
    if 'Gene' in df.columns:
        df['Gene'] = df['Gene'].str.upper()

    df['Mask'] = mask_label
    return df


####################################################
####### DATA INITIALISATION ########################
####################################################

# Load Olink target proteins
proteins_df = pd.read_csv(f'{DATA_DIR}/Olink Proteins.csv', index_col=False)
proteins_upper = set(proteins_df['PROTEIN'].str.upper())
proteins_lower = set(proteins_df['PROTEIN'].str.lower())

# Load UKB variants
clean_mane = pd.read_parquet(
    '/slade/home/mr935/scripts/biobank_work/parquets/clean_mane.parquet',
    engine='pyarrow',
)

####################################################
####### MISSENSE SPLICEOGENIC VARIANTS #############
####################################################
candidate_vars = clean_mane[
    (clean_mane['Max_SpliceAI_Score'] >= 0.8)
    & (clean_mane['REVEL'] < 0.5)
    & (clean_mane['am_pathogenicity'] < 0.34)
]
candidate_vars.to_csv(
    f'{BIOBANK_DIR}/candidate_vars.csv', index=False, na_rep='-'
)

olink_vars = candidate_vars[candidate_vars['SYMBOL'].isin(proteins_upper)].copy()
olink_vars.to_csv(f'{BIOBANK_DIR}/olink_candidate_vars.csv', index=False)

format_and_export_mask(
    olink_vars,
    gene_col='SYMBOL',
    var_col='variantID',
    output_path=f'{BIOBANK_DIR}/check_olink_variants.tsv',
)

# Load extracted Olink data
olink_data = clean_olink_output(
    f'{BIOBANK_DIR}/olink_var_data.tsv', 'p.Splice Missense'
)
olink_data.to_csv(f'{BIOBANK_DIR}/my_olink_data.tsv', sep='\t', index=False)

####################################################
####### MISSENSE DAMAGING VARIANTS #################
####################################################
missense_damaging = clean_mane[
    (clean_mane['REVEL'] >= 0.932)
    & (clean_mane['Max_SpliceAI_Score'] <= 0.01)
    & (clean_mane['am_pathogenicity'] >= 0.564)
    & (clean_mane['SYMBOL'].str.lower().isin(proteins_lower))
].copy()

format_and_export_mask(
    missense_damaging,
    gene_col='SYMBOL',
    var_col='variantID',
    output_path=f'{BIOBANK_DIR}/missense_mask.tsv',
)

# Load extracted Olink data
olink_missense = clean_olink_output(
    f'{OLINK_DIR}/olink_missense_data.tsv', 'p.Missense Damaging'
)
olink_missense.to_csv(
    f'{BIOBANK_DIR}/missense_mask.tsv', sep='\t', index=False
)

####################################################
####### MISSENSE CONTROL VARIANTS ##################
####################################################
missense_benign = clean_mane[
    clean_mane['Consequence'].str.contains('missense', na=False)
    & (clean_mane['Max_SpliceAI_Score'] <= 0.01)
    & (clean_mane['REVEL'] < 0.5)
    & (clean_mane['am_pathogenicity'] < 0.34)
    & (clean_mane['CADD_PHRED'] <= 0.15)
    & (clean_mane['SYMBOL'].str.lower().isin(proteins_lower))
].copy()

format_and_export_mask(
    missense_benign,
    gene_col='SYMBOL',
    var_col='variantID',
    output_path=f'{BIOBANK_DIR}/missense_benign_mask.tsv',
)

# Load extracted Olink data
olink_missense_benign = clean_olink_output(
    f'{OLINK_DIR}/olink_missense_benign_data.tsv', 'Missense Control'
)
olink_missense_benign.to_csv(
    f'{BIOBANK_DIR}/missense_benign_olink_mask.tsv', sep='\t', index=False
)

####################################################
####### SYNONYMOUS CONTROL VARIANTS ################
####################################################
synonymous_benign = clean_mane[
    clean_mane['Consequence'].str.contains('synonymous', na=False)
    & (clean_mane['CADD_PHRED'] <= 0.15)
    & (clean_mane['Max_SpliceAI_Score'] <= 0.01)
    & (clean_mane['SYMBOL'].str.lower().isin(proteins_lower))
].copy()

format_and_export_mask(
    synonymous_benign,
    gene_col='SYMBOL',
    var_col='variantID',
    output_path=f'{BIOBANK_DIR}/synonymous_benign_prots.tsv',
)

# Load extracted Olink data
olink_syn_ben = clean_olink_output(
    f'{BIOBANK_DIR}/olink_synonymous_benign_data.tsv', 'Synonymous Control'
)
olink_syn_ben.to_csv(
    f'{BIOBANK_DIR}/olink_synonymous_benign_mask.tsv', sep='\t', index=False
)

####################################################
####### SYNONYMOUS SPLICEOGENIC VARIANTS ###########
####################################################
synonymous = clean_mane[
    clean_mane['Consequence'].str.contains('synonymous', na=False)
]
splice_syn = synonymous[
    (synonymous['Max_SpliceAI_Score'] >= 0.8)
    & (synonymous['SYMBOL'].str.lower().isin(proteins_lower))
].copy()

format_and_export_mask(
    splice_syn,
    gene_col='SYMBOL',
    var_col='variantID',
    output_path=f'{BIOBANK_DIR}/splice_synonymous_prots.tsv',
)

# Load extracted Olink data
splice_syn_mask = clean_olink_output(
    f'{BIOBANK_DIR}/olink_splice_synonymous_prots_data.tsv',
    'p.Splice Synonymous',
)

####################################################
####### LOSS-OF-FUNCTION VARIANTS ##################
####################################################
# Load VEP output using comment='#' to handle variable metadata lines cleanly
lofs = pd.read_csv(
    f'{BIOBANK_DIR}/LoF_ukb_variants.vep.txt', sep='\t', comment='#'
)
cadd_lofs = lofs[
    (lofs['CADD_PHRED'].astype(float) >= 28.1)
    & (lofs['SYMBOL'].str.lower().isin(proteins_lower))
].copy()

# Save input mask for bash extraction
cadd_lofs.to_csv(f'{BIOBANK_DIR}/lof_variants.tsv', sep='\t', index=False)

# Load extracted Olink data & filter indels
cadd_lofs_mask = clean_olink_output(
    f'{BIOBANK_DIR}/olink_lof_data.tsv', 'p.LoF'
)
cadd_lofs_mask.drop_duplicates(inplace=True)

# Filter out insertions and deletions based on allele lengths
variant_split = cadd_lofs_mask['Variant'].str.split(':', expand=True)
non_indel = (variant_split[2].str.len() <= 1) & (
    variant_split[3].str.len() <= 1
)
cadd_lofs_mask = cadd_lofs_mask[non_indel].reset_index(drop=True)

cadd_lofs_mask.to_csv(
    f'{BIOBANK_DIR}/olink_lof_mask.tsv', sep='\t', index=False
)

####################################################
####### SPLICEVARDB VARIANTS #######################
####################################################
splice_vars = pd.read_csv(
    f'{DATA_DIR}/splicevardb_full_set.tsv', sep='\t', index_col=False
)
splice_vars_filtered = splice_vars[
    splice_vars['gene'].str.lower().isin(proteins_lower)
].copy()

format_and_export_mask(
    splice_vars_filtered,
    gene_col='gene',
    var_col='hg38',
    output_path=f'{BIOBANK_DIR}/splice_var_set.tsv',
)

# Load extracted Olink data
olink_splice = clean_olink_output(
    f'{BIOBANK_DIR}/olink_splice_var_full_data.tsv', 'Validated Splice'
)
olink_splice.to_csv(
    f'{BIOBANK_DIR}/olink_splice_mask.tsv', sep='\t', index=False
)

####################################################
####### COMPILE AND SAVE ###########################
####################################################
compiled_olink_data = pd.concat(
    [
        olink_data,
        olink_missense,
        olink_missense_benign,
        olink_syn_ben,
        splice_syn_mask,
        cadd_lofs_mask,
        olink_splice,
    ],
    ignore_index=True,
)

compiled_olink_data.to_csv(
    f'{BIOBANK_DIR}/olink_files/final_olink_data.csv', index=False
)
