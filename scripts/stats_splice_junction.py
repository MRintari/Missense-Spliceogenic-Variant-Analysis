import argparse
import sys
import os
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.contingency_tables import StratifiedTable


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run Intra- and Inter-Sample splicing statistical tests on read-level CSV metadata."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to the input CSV file containing read-level splicing data."
    )
    return parser.parse_args()


def run_statistical_tests(csv_file):
    if not os.path.exists(csv_file):
        sys.exit(f"Error: CSV file '{csv_file}' does not exist.")

    df = pd.read_csv(csv_file)
    print(f"Loaded {len(df)} read entries from '{csv_file}'.\n")

    # ==============================================================================
    # TEST 1: INTRA-SAMPLE TEST (ALT vs REF inside Carriers)
    # ==============================================================================
    carrier_df = df[df['Group'] == 'Carrier']
    
    # Filter for reads where an allele was actually assigned (ALT or REF)
    carrier_allele_df = carrier_df[carrier_df['Allele'].isin(['ALT', 'REF'])]

    print("=" * 60)
    print("TEST 1: INTRA-SAMPLE ANALYSIS (Carriers ALT vs REF)")
    print("=" * 60)

    if len(carrier_allele_df) == 0:
        print("Warning: No carrier reads with assigned ALT/REF alleles found.")
    else:
        intra_ct = pd.crosstab(carrier_allele_df['Allele'], carrier_allele_df['Splicing_Status'])
        print("\nContingency Table (Carriers Only):")
        print(intra_ct)

        # Ensure correct orientation for Fisher's test
        odds_ratio, p_val_intra = stats.fisher_exact(intra_ct)
        print(f"\nFisher's Exact Test Results:")
        print(f"  Odds Ratio: {odds_ratio:.4f}")
        print(f"  p-value:    {p_val_intra:.4e}")

    # ==============================================================================
    # TEST 2A: STRATIFIED CMH TEST (Carriers vs Controls, Stratified by Date)
    # ==============================================================================
    print("\n" + "=" * 60)
    print("TEST 2A: INTER-SAMPLE ANALYSIS (CMH Test Stratified by Date)")
    print("=" * 60)

    dates = df['Seq_Date'].unique()
    tables = []

    for d in dates:
        sub_df = df[df['Seq_Date'] == d]
        ct = pd.crosstab(sub_df['Group'], sub_df['Splicing_Status'])
        # Reindex to ensure 2x2 shape even if a group has 0 counts in a specific batch
        ct = ct.reindex(index=['Carrier', 'Control'], columns=['Retained', 'Spliced'], fill_value=0)
        print(ct)
        tables.append(ct.values)

    try:
        cmh = StratifiedTable(tables)
        print(f"\nUnique Sequencing Dates Identified: {list(dates)}")
        print(f"Pooled Odds Ratio: {cmh.oddsratio_pooled:.4f}")
        print(f"CMH Test p-value:  {cmh.test_null_odds().pvalue:.4e}")
    except Exception as e:
        print(f"Could not calculate CMH test: {e}")

    # ==============================================================================
    # TEST 2B: MIXED-EFFECTS LOGISTIC REGRESSION (GLMM)
    # ==============================================================================
    print("\n" + "=" * 60)
    print("TEST 2B: INTER-SAMPLE ANALYSIS (Mixed-Effects Logistic Regression)")
    print("=" * 60)

    # Set Control as reference group
    df['Group'] = pd.Categorical(df['Group'], categories=['Control', 'Carrier'])

    try:
        # Fits: Retention ~ Group + Sequencing_Date with Random Effect on Sample_ID
        glmm = smf.mixedlm("Retained_Binary ~ C(Group) + C(Seq_Date)", df, groups=df["Sample_ID"])
        glmm_result = glmm.fit()
        print("\nModel Summary:")
        print(glmm_result.summary())
    except Exception as e:
        print(f"Could not fit Mixed-Effects GLMM model: {e}")
        print("Fallback: Fitting standard Logistic Regression (GLM) without random effects...")
        glm = smf.glm("Retained_Binary ~ C(Group) + C(Seq_Date)", df, family=sm.families.Binomial()).fit()
        print(glm.summary())

    print("\n" + "=" * 60)


def main():
    args = parse_arguments()
    run_statistical_tests(args.input)


if __name__ == "__main__":
    main()
