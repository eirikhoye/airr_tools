import pandas as pd
import sys
import os

"""
Script to concatenate rearrangement files, 

Usage: python combine_cdr3_files.py <rearrangement_file_1.tsv> <rearrangement_file_2.tsv> <rearrangement_file_[...].tsv>
"""

path_list = sys.argv[1:-1]

is_first = True
for file in path_list:
    if is_first:
        df = pd.read_csv(file, sep='\t')
        is_first = False
    else:
        df = pd.concat([df, pd.read_csv(file, sep='\t')])

df.to_csv(sys.stdout, index=False, sep='\t')
