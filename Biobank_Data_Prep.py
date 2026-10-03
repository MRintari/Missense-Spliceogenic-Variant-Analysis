import pandas as pd

# Load re-annotated Biobank data
mane = pd.read_csv('/slade/home/mr935/data/biobank_data/samples/updated/splice_mane_filtered_ukb24310.vep.tsv', sep="\t", na_values=['-'], index_col=False)
cx = pd.read_csv('/slade/home/mr935/data/biobank_data/samples/updated/splice_mane_filtered_ukb24310_cX.vep.tsv', sep="\t", na_values=['-'], index_col=False)
cy = pd.read_csv('/slade/home/mr935/data/biobank_data/samples/updated/splice_mane_filtered_ukb24310_cY.vep.tsv', sep="\t", na_values=['-'], index_col=False)

cx = cx[['#Uploaded_variation', 'SpliceAI_Pred']]
cy = cy[['#Uploaded_variation', 'SpliceAI_Pred']]

# Reformat ChrX,Y to match autosomal data
fix_mane = pd.merge(mane, cx, on='#Uploaded_variation', how='left', suffixes=('_original', '_correct'))
fix_mane.loc[fix_mane['variantID'].str.contains('chrX'), 'SpliceAI_pred'] = fix_mane.loc[fix_mane['variantID'].str.contains('chrX'), 'SpliceAI_Pred']
fix_mane.drop('SpliceAI_Pred', axis=1, inplace=True)

fix_mane = pd.merge(fix_mane, cy, on='#Uploaded_variation', how='left', suffixes=('_original', '_correct'))
fix_mane.loc[fix_mane['variantID'].str.contains('chrY'), 'SpliceAI_pred'] = fix_mane.loc[fix_mane['variantID'].str.contains('chrY'), 'SpliceAI_Pred']
fix_mane.drop('SpliceAI_Pred', axis=1, inplace=True)

clean_mane = fix_mane.copy()

# Clean annotations to reflect the MANE SELECT transcript
def clean_spliceAI(row):
    splice = row['SpliceAI_pred']
    feature = row['Feature']
    if pd.notna(splice) and splice:
        results = splice.split(',')
        for res in results:
            scores = res.split('|')
            transcript = scores[1].split('.')[0]
            if feature == transcript:
                return res
    return splice

clean_mane['SpliceAI_pred'] = clean_mane.apply(clean_spliceAI, axis=1)

# Filter out genes missing valid SpliceAI scores
invalid_mask = clean_mane['SpliceAI_pred'].astype(str).str.split('|').str.len() < 3
genes_to_filter = clean_mane[invalid_mask]['SYMBOL'].unique()

clean_mane = clean_mane[~(clean_mane['SYMBOL'].isin(genes_to_filter))]

clean_mane.to_parquet('clean_mane.parquet')

