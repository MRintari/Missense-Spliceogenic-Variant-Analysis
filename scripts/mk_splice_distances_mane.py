# Import Modules
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pyranges as pr

# Import annotations
ann = pr.read_gtf('/slade/home/mr935/SpliceAI-lookup/annotations/gencode.v48.basic.annotation.gtf.gz')
ann_df = ann.df

# Select exons
exons = ann[(ann.Feature == 'exon')]
exons = exons[['gene_id', 'gene_name', 'transcript_id']]

# MANE Select Transcripts
mane = pd.read_csv('/slade/home/mr935/SpliceAI-lookup/annotations/MANE.GRCh38.v1.4.transcripts_by_gene.tsv', sep="\t", index_col=False)
mane_exons = exons[exons.transcript_id.isin(mane['MANE_Select_Ensembl_id'])]
mane_exons_df = mane_exons.df

# Identify which transcripts need to be extended and at what ends
def mark_splice_sites(df):
    # For each exon in MANE Select transcripts
    def process_transcript(transcript):
        transcript = transcript.reset_index(drop=True)
        n = len(transcript)

        # Initialize flags to create donor and acceptor sites
        transcript['make_acceptor'] = True
        transcript['make_donor'] = True

        # First exon: no acceptor site
        transcript.loc[0, 'make_acceptor'] = False
        # Last exon: no donor site
        transcript.loc[n-1, 'make_donor'] = False

        # If single exon transcript, no donor or acceptor
        if n == 1:
            transcript['make_acceptor'] = False
            transcript['make_donor'] = False

        return transcript

    return df.groupby('transcript_id').apply(process_transcript).reset_index(drop=True)

exons_marked_df = mark_splice_sites(mane_exons_df)

# Split by Splice Site Type
acceptors_df = exons_marked_df[exons_marked_df['make_acceptor']].copy()
donors_df = exons_marked_df[exons_marked_df['make_donor']].copy()

# Clean to enable PyRanges conversion
acceptors_df['Strand'] = acceptors_df['Strand'].cat.remove_unused_categories()
donors_df['Strand'] = donors_df['Strand'].cat.remove_unused_categories()
exons_marked_df['Strand'] = exons_marked_df['Strand'].cat.remove_unused_categories()

exons_marked = pr.PyRanges(exons_marked_df)
acceptors = pr.PyRanges(acceptors_df)
donors = pr.PyRanges(donors_df)

# Get exon's acceptor/donor sites
acceptor_ints = acceptors.extend({"5": 1, "3": 0}).subtract(acceptors).extend({"5": 0, "3": 1}).subsequence(-1)
donor_ints = donors.extend({"5": 0, "3": 1}).subtract(donors)

acceptor_ints = acceptor_ints.assign("SpliceSite", lambda df: pd.Series(["Acceptor"] * len(df)))
donor_ints = donor_ints.assign("SpliceSite", lambda df: pd.Series(["Donor"] * len(df)))

acceptor_ints_df = acceptor_ints.df
donor_ints_df = donor_ints.df

splice_sites_df = pd.concat([acceptor_ints_df, donor_ints_df])
splice_sites_df.drop(['make_acceptor', 'make_donor'], inplace=True, axis=1)

splice_sites = pr.PyRanges(splice_sites_df)

# Save GTF output
splice_sites.to_gtf('/slade/home/mr935/scripts/Splice_Distance/data/mane_splice_sites.gencode.v48.basic.annotation.gtf.gz')
