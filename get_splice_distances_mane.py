# Import modules
import pyranges as pr
import pandas as pd

# Load Data
clean_mane = pd.read_parquet('clean_mane.parquet')
mane_splice_sites = pr.read_gtf('/slade/home/mr935/scripts/Splice_Distance/data/mane_splice_sites.gencode.v48.basic.annotation.gtf.gz')

# Generate Donor and Acceptor Annotation Files
donors = mane_splice_sites[mane_splice_sites.SpliceSite == 'Donor']
acceptors = mane_splice_sites[mane_splice_sites.SpliceSite == 'Acceptor']

# Generate necessary PyRanges columns
clean_mane['Chromosome'] = clean_mane['variantID'].str.split("-").str[0]
clean_mane['Position'] = clean_mane['variantID'].str.split("-").str[1].astype(int)
clean_mane['Strand'] = np.where(clean_mane['STRAND'] == 1, '+', '-')
clean_mane['Start'] = clean_mane['Position'] - 1
clean_mane['End'] = clean_mane['Position'] 

# Identify nearest Donor and Acceptor Sites per variant
clean_mane_pr = pr.PyRanges(clean_mane)
nearest_donors = clean_mane_pr.nearest(donors)
nearest_acceptors = clean_mane_pr.nearest(acceptors)

# Identify if Acceptor or Donor Site is closer via merge and .where()
donors_df = nearest_donors.df
acceptors_df = nearest_acceptors.df

clean_mane = clean_mane.set_index('variantID')
donors_df = donors_df.set_index('variantID')
acceptors_df = acceptors_df.set_index('variantID')

merge_cols = ['Distance', 'Start_b', 'End_b']

clean_mane = clean_mane.merge(
    donors_df[merge_cols],
    left_index=True,  # Use the index variantID
    right_index=True, 
    how='left',
    suffixes=('', '_Donor')
)

clean_mane = clean_mane.merge(
    acceptors_df[merge_cols],
    left_index=True,
    right_index=True,
    how='left',
    suffixes=('_Donor', '_Acceptor')
)

clean_mane['Nearest_Site'] = np.where(clean_mane['Distance_Donor'] <
                                      clean_mane['Distance_Acceptor'],
                                      'Donor', 'Acceptor')

# Drop PyRange Columns
clean_mane.drop(['Distance_Donor', 'Start_b_Donor', 'End_b_Donor',
                 'Distance_Acceptor', 'Start_b_Acceptor', 'End_b_Acceptor'],
                 axis=1, inplace=True)


# Reset variant intervals for accurate distance generation
is_donor = clean_mane['Nearest_Site'] == 'Donor'

# Acceptor: [Position - 1, Position] (Already in 'Start' and 'End')
# Donor: [Position, Position + 1]
clean_mane['Start'] = np.where(is_donor, clean_mane['Position'], clean_mane['Start'])
clean_mane['End'] = np.where(is_donor, clean_mane['Position'] + 1, clean_mane['End'])

# Generate accurate distances for each variant
clean_mane_pr = pr.PyRanges(clean_mane)
nearest = clean_mane_pr.nearest(mane_splice_sites)

# Merge accurate distances to mane DF
nearest_df = nearest.df
clean_mane = clean_mane.reset_index()
clean_mane = clean_mane.set_index('#Uploaded_variation')
nearest_df = nearest_df.set_index('#Uploaded_variation')


clean_mane = clean_mane.merge(
    nearest_df[merge_cols],
    left_index=True,
    right_index=True,
    how='left'
)

# Clean DF
clean_mane = clean_mane.reset_index()
clean_mane = clean_mane.rename(columns={'Distance': 'Variant_Canonical_Splice_Distance',
                                        'Start_b': 'Canonical_Splice_Int_Start',
                                        'End_b': 'Canonical_Splice_Int_End'})

clean_mane.to_parquet('clean_mane.parquet')
