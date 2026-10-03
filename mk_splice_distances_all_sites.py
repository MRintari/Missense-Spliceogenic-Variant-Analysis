# Import Modules
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pyranges as pr

# Note: Unlike mk_splice_distances_mane.ipynb, we cannot use the .subtract() method
# on all transcripts. Transcripts with splice sites within an alternate transcript's 
# exon are not found as they are subtracted out when using .subtract().


def mark_splice_sites(df):
    """Identify which transcripts need to be extended and at what ends."""
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


# ==========================================
# 1. Load & Process GENCODE Exons
# ==========================================
gencode = pr.read_gtf('/slade/home/mr935/SpliceAI-lookup/annotations/gencode.v48.basic.annotation.gtf.gz')

exons = gencode[gencode.Feature == 'exon']
exons = exons[['gene_id', 'gene_name', 'transcript_id']]
exons_df = exons.df

exons_marked_df = mark_splice_sites(exons_df)

# Split by Splice Site Type and Clean
acceptors_df = exons_marked_df[exons_marked_df['make_acceptor']].copy()
donors_df = exons_marked_df[exons_marked_df['make_donor']].copy()

acceptors_df['Strand'] = acceptors_df['Strand'].cat.remove_unused_categories()
donors_df['Strand'] = donors_df['Strand'].cat.remove_unused_categories()

# Split DFs by Strand to establish interval rules
acc_pos = acceptors_df[acceptors_df['Strand'] == '+'].copy().reset_index(drop=True)
acc_neg = acceptors_df[acceptors_df['Strand'] == '-'].copy().reset_index(drop=True)
don_pos = donors_df[donors_df['Strand'] == '+'].copy().reset_index(drop=True)
don_neg = donors_df[donors_df['Strand'] == '-'].copy().reset_index(drop=True)

# Identify Splice Sites
acc_pos['Start_b'] = acc_pos['Start']
acc_pos['End_b'] = acc_pos['Start'] + 1

acc_neg['Start_b'] = acc_neg['End']
acc_neg['End_b'] = acc_neg['End'] + 1

don_pos['Start_b'] = don_pos['End']
don_pos['End_b'] = don_pos['End'] + 1

don_neg['Start_b'] = don_neg['Start']
don_neg['End_b'] = don_neg['Start'] + 1

# Compile and Clean DFs
donors = pd.concat([don_neg, don_pos])
donors['SpliceSite'] = pd.Series(['Donor'] * len(donors))

acceptors = pd.concat([acc_neg, acc_pos])
acceptors['SpliceSite'] = pd.Series(['Acceptor'] * len(acceptors))

for df in [acceptors, donors]:
    df.drop(columns=['Start', 'End'], axis=1, inplace=True)
    df.rename(columns={'Start_b': 'Start', 'End_b': 'End'}, inplace=True)

gencode_splices = pd.concat([donors, acceptors])
gencode_splices.drop(columns=['make_acceptor', 'make_donor'], axis=1, inplace=True)
gencode_splices.drop('gene_id', axis=1, inplace=True)
gencode_splices.rename(columns={'gene_name': 'gene'}, inplace=True)


# ==========================================
# 2. Load & Process RefSeq Exon Annotations
# ==========================================
refseq = pr.read_gff3('/slade/home/mr935/data/refseq_ncbi.GRCh38.p14.annotation.features.gff.gz')

exons = refseq[refseq.Feature == 'exon']
exons = exons[['gene', 'transcript_id']]
exons_df = exons.df

exons_marked_df = mark_splice_sites(exons_df)

# Split by Splice Site Type and Clean
acceptors_df = exons_marked_df[exons_marked_df['make_acceptor']].copy()
donors_df = exons_marked_df[exons_marked_df['make_donor']].copy()

acceptors_df['Strand'] = acceptors_df['Strand'].cat.remove_unused_categories()
donors_df['Strand'] = donors_df['Strand'].cat.remove_unused_categories()

# Split DFs by Strand to establish interval rules
acc_pos = acceptors_df[acceptors_df['Strand'] == '+'].copy().reset_index(drop=True)
acc_neg = acceptors_df[acceptors_df['Strand'] == '-'].copy().reset_index(drop=True)
don_pos = donors_df[donors_df['Strand'] == '+'].copy().reset_index(drop=True)
don_neg = donors_df[donors_df['Strand'] == '-'].copy().reset_index(drop=True)

# Identify Splice Sites
acc_pos['Start_b'] = acc_pos['Start']
acc_pos['End_b'] = acc_pos['Start'] + 1

acc_neg['Start_b'] = acc_neg['End']
acc_neg['End_b'] = acc_neg['End'] + 1

don_pos['Start_b'] = don_pos['End']
don_pos['End_b'] = don_pos['End'] + 1

don_neg['Start_b'] = don_neg['Start']
don_neg['End_b'] = don_neg['Start'] + 1

# Compile and Clean DFs
donors = pd.concat([don_neg, don_pos])
donors['SpliceSite'] = pd.Series(['Donor'] * len(donors))

acceptors = pd.concat([acc_neg, acc_pos])
acceptors['SpliceSite'] = pd.Series(['Acceptor'] * len(acceptors))

for df in [acceptors, donors]:
    df.drop(columns=['Start', 'End'], axis=1, inplace=True)
    df.rename(columns={'Start_b': 'Start', 'End_b': 'End'}, inplace=True)

refseq_splices = pd.concat([donors, acceptors])
refseq_splices.drop(columns=['make_acceptor', 'make_donor'], axis=1, inplace=True)

# Standardise RefSeq Chrom Names
refseq_names = {
    'NC_000001': 'chr1',  'NC_000002': 'chr2',  'NC_000003': 'chr3',
    'NC_000004': 'chr4',  'NC_000005': 'chr5',  'NC_000006': 'chr6',
    'NC_000007': 'chr7',  'NC_000008': 'chr8',  'NC_000009': 'chr9',
    'NC_000010': 'chr10', 'NC_000011': 'chr11', 'NC_000012': 'chr12',
    'NC_000013': 'chr13', 'NC_000014': 'chr14', 'NC_000015': 'chr15',
    'NC_000016': 'chr16', 'NC_000017': 'chr17', 'NC_000018': 'chr18',
    'NC_000019': 'chr19', 'NC_000020': 'chr20', 'NC_000021': 'chr21',
    'NC_000022': 'chr22', 'NC_000023': 'chrX',  'NC_000024': 'chrY'
}

def rename_refseq_chroms(contig):
    if contig.startswith('NC_'):
        chrom = contig.split('.')[0]
        return refseq_names.get(chrom, contig)
    return contig

refseq_splices['Chromosome'] = refseq_splices['Chromosome'].apply(rename_refseq_chroms)


# ==========================================
# 3. Compile Final Splice Site Annotations
# ==========================================
all_splices = pd.concat([gencode_splices, refseq_splices])

# Export output GFF3 file
pr.PyRanges(all_splices).to_gff3('/slade/home/mr935/data/all_annotated_splice_sites.v48.GRCh38.annotation.gff.gz')