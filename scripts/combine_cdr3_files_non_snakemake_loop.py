import pandas as pd
import sys
import os


path = "data/gliph_data/"
path_list = os.listdir(path)


is_first = True
for file in path_list:
    if is_first:
        df = pd.read_csv(path+file, sep='\t')
        is_first = False
    else:
        df = pd.concat([df, pd.read_csv(path+file, sep='\t')])

df.to_csv("data/combined_cdr3_file.txt", index=False, sep='\t')