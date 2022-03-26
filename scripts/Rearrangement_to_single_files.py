import pandas as pd
import sys

path = sys.argv[1]
in_file = sys.argv[2]
df = pd.read_csv(path+in_file, sep='\t')

file_list = df['sample_name'].unique()
for file in file_list:
    print(file)
    file_name = file+'.tsv'
    new_df = df[df['sample_name'] == file]
    new_df.to_csv(path+file_name, sep='\t')