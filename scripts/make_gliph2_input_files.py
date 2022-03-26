import pandas as pd
import sys
import os

# path to ImmunSEQ rearrangement file
path = sys.argv[1]

# Read pandas dataframe
df = pd.read_csv(path, sep='\t')
df = df[['amino_acid', 'v_resolved', 'j_resolved', 'sample_name', 'templates']]
df.columns = ['CDR3b', 'TRBV', 'TRBJ', 'subject', 'count']

# Reformat V_rearrangement to GLIPH2 format
df['TRBV'] = [i.split('*')[0] for i in df['TRBV']]
df['TRBV'] = [i.replace('TCRBV', 'TRBV') for i in df['TRBV']]
df['TRBV'] = [i.replace('TRBV0', 'TRBV') for i in df['TRBV']]
df['TRBV'] = [i.replace('-0', '-') for i in df['TRBV']]

# Reformat J_rearrangement to GLIPH2 format
df['TRBJ'] = [i.split('*')[0] for i in df['TRBJ']]
df['TRBJ'] = [i.replace('TCRBJ', 'TRBJ') for i in df['TRBJ']]
df['TRBJ'] = [i.replace('TRBJ0', 'TRBJ') for i in df['TRBJ']]
df['TRBJ'] = [i.replace('-0', '-') for i in df['TRBJ']]

# Write tab separated to stdout
df.to_csv(sys.stdout, index=False, sep='\t')