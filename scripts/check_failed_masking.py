# Import modules
import pyranges as pr
import pandas as pd

# Load MANE Select Splice Junctions
mane_splice_sites = pr.read_gtf('/slade/home/mr935/scripts/Splice_Distance/data/mane_splice_sites.gencode.v48.basic.annotation.gtf.gz')
mane_donors = mane_splice_sites[mane_splice_sites.SpliceSite == 'Donor']
mane_acceptors = mane_splice_sites[mane_splice_sites.SpliceSite == 'Acceptor']

# Load LMSVs
candidate_vars = pd.read_csv('/slade/home/mr935/data/biobank_data/candidate_vars.csv', index_col=False)

# Get intervals for start,end of splice sites
# Ensures format of site start,ends match for both donor and acceptor sites 
def get_ints(row):
    effect = row['Max_SpliceAI_Column']
    if 'DL' in effect or 'DG' in effect:
        return [row['Splice_Pos'], row['Splice_Pos'] + 1]
    else:
        return [row['Splice_Pos'] - 1, row['Splice_Pos']]
    
# Prepare candidate_vars data for PyRanges conversion
candidate_vars['Chromosome'] = candidate_vars['variantID'].str.split("-").str[0]
candidate_vars['Splice_Pos'] = pd.to_numeric(candidate_vars['variantID'].str.split("-").str[1]) + candidate_vars['Max_SpliceAI_Position']
candidate_vars['Splice_Pos'] = candidate_vars['Splice_Pos'].astype(int)
candidate_vars['Ints'] = candidate_vars.apply(get_ints, axis=1)
candidate_vars['Start'] = candidate_vars['Ints'].str[0]
candidate_vars['End'] = candidate_vars['Ints'].str[1]
candidate_vars.drop('Ints', axis=1, inplace=True)

# Function to determine if masking functioned correctly on 
# variant of interest. Assesses if the delta position lands
# on a MANE Select splice junction, i.e. variant predicted
# to increase strength of transcript's own splice junction
def check_masking(row):
    # Convert row to PyRange object
    row = pd.DataFrame([row])
    pyrange = pr.PyRanges(row)

    # Determine Distance dependent on Splice Effect
    effect = list(pyrange.Max_SpliceAI_Column)[0]
    match effect:
        case 'SpliceAI_DS_DL':
            return 'Loss Event'
        case 'SpliceAI_DS_AL':
            return 'Loss Event'
        case 'SpliceAI_DS_DG':
            distance = pyrange.nearest(mane_donors)
            distance = list(distance.Distance)[0]
            return True if distance == 0 else False
        case 'SpliceAI_DS_AG':
            distance = pyrange.nearest(mane_donors)
            distance = list(distance.Distance)[0]
            return True if distance == 0 else False
        case _:
            return None
        
candidate_vars['Bad_Mask?'] = candidate_vars.apply(check_masking, axis=1)
len(candidate_vars[candidate_vars['Bad_Mask?'] == True]) # Returns 0, i.e no bad masking present
