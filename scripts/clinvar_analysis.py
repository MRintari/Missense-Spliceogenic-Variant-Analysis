# Import modules
import pandas as pd

# Load ClinVar Data
clean_benign_missense = pd.read_parquet('/slade/home/mr935/data/clinvar_annots/clean_missense_annots/clean_benign_missense.parquet')
clean_lb_missense = pd.read_parquet('/slade/home/mr935/data/clinvar_annots/clean_missense_annots/clean_lb_missense.parquet')
clean_vus_missense = pd.read_parquet('/slade/home/mr935/data/clinvar_annots/clean_missense_annots/clean_vus_missense.parquet')
clean_lp_missense = pd.read_parquet('/slade/home/mr935/data/clinvar_annots/clean_missense_annots/clean_lp_missense.parquet')
clean_pathogenic_missense = pd.read_parquet('/slade/home/mr935/data/clinvar_annots/clean_missense_annots/clean_pathogenic_missense.parquet')


def max_splice(df):
    """Extract Variant's Highest SpliceAI Score"""
    splice_deltas = ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']
    df = df.copy()

    df[splice_deltas] = df[splice_deltas].astype(float)
    df['Max_SpliceAI'] = df[splice_deltas].max(axis=1)

    return df['Max_SpliceAI']

clean_benign_missense['Max_SpliceAI'] = max_splice(clean_benign_missense)
clean_lb_missense['Max_SpliceAI'] = max_splice(clean_lb_missense)
clean_vus_missense['Max_SpliceAI'] = max_splice(clean_vus_missense)
clean_lp_missense['Max_SpliceAI'] = max_splice(clean_lp_missense)
clean_pathogenic_missense['Max_SpliceAI'] = max_splice(clean_pathogenic_missense)

sig_dfs = [clean_benign_missense, clean_lb_missense, clean_vus_missense,
           clean_lp_missense, clean_pathogenic_missense]
clnsigs = ['Benign', 'Likely Benign', 'VUS', 'Likely Pathogenic', 'Pathogenic']

def get_percs(df, clinvar):
    size = len(df)
    two = len(df[df['Max_SpliceAI'] >= 0.2]) / size * 100
    five = len(df[df['Max_SpliceAI'] >= 0.5]) / size * 100
    eight = len(df[df['Max_SpliceAI'] >= 0.8]) / size * 100

    return print(f'Of the {size} {clinvar} variants:\n{two:.2f}%
                 have a SpliceAI score \u2265 0.2\n{five:.2f}% have a SpliceAI score
                 \u2265 0.5\n{eight:.2f}% have a SpliceAI score \u2265 0.8\n')

for sig_df, sig in zip(sig_dfs, clnsigs):
    get_percs(sig_df, sig)
