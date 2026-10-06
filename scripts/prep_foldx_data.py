# Import modules
import pandas as pd

# Load data
foldx = pd.read_csv('/slade/home/mr935/data/2025.02.10_foldx_energy.csv.gz', index_col=False)
olink_data = pd.read_csv('/slade/home/mr935/data/biobank_data/olink_files/final_olink_data.csv',
                         index_col=False)
clean_mane = pd.read_parquet('biobank_work/parquets/clean_mane.parquet', engine='pyarrow')

# Get the uniprot accession IDs for each gene
uniprot = pd.read_csv('/slade/home/mr935/data/uniprot_ids.tsv', sep='\t', index_col=False)

# Format
olink_data['variantID'] = 'chr' + olink_data['Variant'].str.replace(':', '-')

# Identify data to plot for FoldX analysis
protvar_checks = olink_data[olink_data['Mask'].
                            isin(['Spliceogenic Missense'
                                  'Missense Damaging'])].copy().reset_index(drop=True)

protvar_checks['Gene'] = protvar_checks['Gene'].str.upper()
protvar_checks = protvar_checks.drop(['Approved symbol', 'Approved name'], axis=1)

# Gather variant data
protvar_checks = protvar_checks.merge(
    clean_mane[['variantID', 'Protein_position']],
    left_on='variantID',
    right_on='variantID',
    how='left'
)   

# Use variant data to identify residues for merging with Uniprot
protvar_checks['wild_type'] = protvar_checks['Amino_acids'].str.split('/').str[0]
protvar_checks['mutated_type'] = protvar_checks['Amino_acids'].str.split('/').str[1]
protvar_checks = protvar_checks[~protvar_checks['UniProt ID(supplied by UniProt)'].isna()]

protvar_checks = protvar_checks.merge(
    foldx[['uniprot_accession', 'uniprot_position',
             'wild_type', 'mutated_type', 'foldx_ddg']],
    left_on=['UniProt ID(supplied by UniProt)', 'Protein_position',
             'wild_type', 'mutated_type'],
    right_on=['uniprot_accession', 'uniprot_position',
              'wild_type', 'mutated_type'],
    how='left'
)

protvar_checks = protvar_checks[~protvar_checks['uniprot_accession'].isna()]

# Filter for destabilising and non-destabilising variants
protvar_checks['Destabilising'] = protvar_checks['foldx_ddg'] >= 2

protvar_checks.to_csv('/slade/home/mr935/data/biobank_data/olink_spliceogenic_and_damaging_missense_protvar.tsv', sep='\t', index=False)
