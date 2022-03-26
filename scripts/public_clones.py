# Import libraries
import pandas as pd
import numpy as np
import os


#########################################################################################################

# Define functions

#########################################################################################################

def get_public_clones(airr_dict):
    """
    Takes as input a dictionary of Adaptive Biotechnologies
    datasets, then loops through all to identify CDR3_aa
    sequences present in more than one sample. Here only
    sequencial information is considered.

    Returns a list of unique public clones.
    """
    public_clones = []

    # Loop through all datasets
    for sample_a in airr_dict:
        print(sample_a)
        seq_a = airr_dict[sample_a].amino_acid.values
        """
        # Compare CDR3_aa sequences in this sample against all
        # other CDR3_aa sequences in all other samples
        """
        for sample_b in airr_dict:
            if sample_b == sample_a:
                continue
            seq_b = airr_dict[sample_b].amino_acid.values
            """
            get shared CDR3_aa
            """
            pairwise_public = np.intersect1d(seq_a, seq_b, assume_unique=False)
            """
            add only new public clones to the public clones list
            """
            public_clones.extend(list(np.setdiff1d(pairwise_public, public_clones)))

        print(len(public_clones))
    return (public_clones)


def get_public_clones_sharing_level(airr_dict):
    """
    Takes as input a dictionary of Adaptive Biotechnologies
    datasets, then loops through all to identify CDR3_aa
    sequences present in more than one sample. Here only
    sequencial information is considered.

    Returns a dict of unique public clones, with public clones
    as keys and nr samples as values
    """
    public_clones = {}

    # Loop through all datasets
    for sample_a in airr_dict:

        print(sample_a)
        seq_a = airr_dict[sample_a].amino_acid.values

        # Compare CDR3_aa sequences in this sample against all
        # other CDR3_aa sequences in all other samples
        for sample_b in airr_dict:
            if sample_b == sample_a:
                continue
            seq_b = airr_dict[sample_b].amino_acid.values

            # get shared CDR3_aa
            pairwise_public = np.intersect1d(seq_a, seq_b, assume_unique=False)

            # add only new public clones to the public clones list
            # public_clones.extend(list(np.setdiff1d(pairwise_public, public_clones)))
            for p_c in pairwise_public:
                if p_c in public_clones:
                    if sample_a not in public_clones[p_c]:
                        public_clones[p_c].append(sample_a)
                    if sample_b not in public_clones[p_c]:
                        public_clones[p_c].append(sample_b)
                else:
                    public_clones[p_c] = [sample_a, sample_b]

        print(len(public_clones))

    public_clones = pd.DataFrame.from_dict({k: len(v) for k, v in public_clones.items()}, orient='index').reset_index()
    public_clones.columns = ['amino_acid', 'sharing_level']
    public_clones = public_clones.sort_values('sharing_level', ascending=False).reset_index(drop=True)
    return (public_clones)


#########################################################################################################

# Gather metadata

#########################################################################################################

# meta_path = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/metadata/Kit_1_metadata.csv'

# meta = pd.read_csv(meta_path, sep=';')
# meta['COMET_ID_group'] = meta['COMET_ID'] + '_group:_' + meta['group']
# meta.index = meta['sample_name']
# meta_dict = meta['COMET_ID_group'].to_dict()


meta_path1 = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/immunarch_data_376276/metadata.txt'
meta1 = pd.read_csv(meta_path1, sep='\t')
meta1['COMET_ID_group'] = meta1['COMET_ID'] + '_group:_' + meta1['group']
meta1.index = meta1['Sample']

meta_path2 = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/immunarch_data_12511620/metadata.txt'
meta2 = pd.read_csv(meta_path2, sep='\t')
meta2['COMET_ID_group'] = meta2['COMET_ID'] + '_group:_' + meta2['group']
meta2.index = meta2['Sample']

meta = pd.concat([meta1, meta2])
meta_dict = meta['COMET_ID_group'].to_dict()

#########################################################################################################

# Make dictionary of airr datasets

#########################################################################################################

# Make a single dictionary of each adaptive biotechnologies airr dataframes, with sample name (LivMet_#)
# as keys and the dataframe as value pair. Input for get_public_clones function.

airr_dict = {}
file_dir = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/immunarch_data_376276/'

for file in os.listdir(file_dir):
    path = file_dir + file
    # if file == 'RearrangementDetails_06-25-2019_12-33-29_AM.tsv':
    #    continue

    if not file.endswith('.tsv'):
        continue
    sample_name = meta_dict[file.split('.')[-2]]
    airr_dict[sample_name] = pd.read_csv(path, sep='\t')

file_dir = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/immunarch_data_12511620/'

for file in os.listdir(file_dir):
    path = file_dir + file
    # if file == 'RearrangementDetails_06-25-2019_12-33-29_AM.tsv':
    #    continue

    if not file.endswith('.tsv'):
        continue
    sample_name = meta_dict[file.split('.')[-2]]
    airr_dict[sample_name] = pd.read_csv(path, sep='\t')

#########################################################################################################

# Use get_public_clones to identify public clones,
# then write to file

#########################################################################################################

out_file_path = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/public_clones.txt'
out_file_path2 = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/public_clones_sharing_level.csv'

remove_samples = [
    'COMET-0011M1_2_group:_no_chemo',
    'COMET-0035M1_2_group:_no_chemo',
    'COMET-0059M1_2_group:_no_chemo',
    'COMET-0032M2_group:_no_chemo',
    'COMET-0026M2_group:_no_chemo',
    'COMET-0030M2_group:_short_chemo',
    'COMET-0062M2_group:_short_chemo',
    'COMET-0041M1_group:_na'
]

airr_dict_sub = {key: airr_dict[key] for key in airr_dict.keys() if key not in remove_samples}

# public_clones = get_public_clones(airr_dict_sub)
# out_file = open(out_file_path, 'w')
# for line in public_clones:
#    out_file.write(line + '\n')

public_clones_sharing_level = get_public_clones_sharing_level(airr_dict_sub)
public_clones_sharing_level.to_csv(out_file_path2, index=False)

out_dir_public = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/rearrangement_files_with_sharing_public/'
for k, v in airr_dict_sub.items():
    public_df = pd.merge(left=v, right=public_clones_sharing_level,
                         left_on='amino_acid', right_on='amino_acid',
                         how='inner')
    public_df = public_df.sort_values(by='templates', ascending=False)
    public_df.to_csv(str(out_dir_public + k + '_public.tsv'), sep='\t')

out_dir_private = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/rearrangement_files_with_sharing_private/'
for k, v in airr_dict_sub.items():
    private_df = v[~v['amino_acid'].isin(public_clones_sharing_level['amino_acid'].values)]
    private_df['sharing_level'] = 1
    private_df = private_df.sort_values(by='templates', ascending=False)
    private_df.to_csv(str(out_dir_private + k + '_private.tsv'), sep='\t')

out_dir_all = '/Users/eirikhoy/Dropbox/projects/airr_comet/data/rearrangement_files_with_sharing_all/'
for k, v in airr_dict_sub.items():
    all_df = pd.merge(left=v, right=public_clones_sharing_level,
                      left_on='amino_acid', right_on='amino_acid',
                      how='left').fillna({'sharing_level': 1})
    all_df = all_df.sort_values(by='templates', ascending=False)
    all_df.to_csv(str(out_dir_all + k + '_all.tsv'), sep='\t')
























































