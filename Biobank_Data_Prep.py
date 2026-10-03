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

# Extract SpliceAI Scores from VEP annotations
def extract_scores(row):
    spliceai = row['SpliceAI_pred']
    spliceai = spliceai.split('|')
    return {
        'SpliceAI_DS_AG': float(spliceai[2]),
        'SpliceAI_DS_AL': float(spliceai[3]),
        'SpliceAI_DS_DG': float(spliceai[4]),
        'SpliceAI_DS_DL': float(spliceai[5]),
        'SpliceAI_DP_AG': int(spliceai[6]),
        'SpliceAI_DP_AL': int(spliceai[7]),
        'SpliceAI_DP_DG': int(spliceai[8]),
        'SpliceAI_DP_DL': int(spliceai[9])
    }

tmp = clean_mane.apply(extract_scores, axis=1)
clean_mane = clean_mane.join(tmp.apply(pd.Series))
clean_mane.drop('SpliceAI_pred', axis=1, inplace=True)

# Find highest splice score
def high_splice(row):
    scores = ['SpliceAI_DS_AG', 'SpliceAI_DS_AL', 'SpliceAI_DS_DG', 'SpliceAI_DS_DL']
    
    return max(list(row[scores]))
clean_mane['Max_SpliceAI_Score'] = clean_mane.apply(high_splice, axis=1)

# Find respective column for high splice score
scores = ['SpliceAI_DS_AG', 'SpliceAI_DS_AL', 'SpliceAI_DS_DG', 'SpliceAI_DS_DL']
score_match = clean_mane[scores] == clean_mane['Max_SpliceAI_Score'].values[:, None]
high_column = score_match.idxmax(axis=1)
clean_mane['Max_SpliceAI_Column'] = np.where(score_match.any(axis=1), high_column, np.nan)

# Find position affected by high splice score
def high_position(row):
    score = row['Max_SpliceAI_Column']
    position_column = score.replace('DS', 'DP')
    position = row[position_column]

    return position

clean_mane['Max_SpliceAI_Position'] = clean_mane.apply(high_position, axis=1)

clean_mane.to_parquet('clean_mane.parquet')
