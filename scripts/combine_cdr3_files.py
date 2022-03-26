import pandas as pd
import sys
import os

path_list = sys.argv[1:-1]

is_first = True
for file in path_list:
    if is_first:
        df = pd.read_csv(file, sep='\t')
        is_first = False
    else:
        df = pd.concat([df, pd.read_csv(file, sep='\t')])

df.to_csv(sys.stdout, index=False, sep='\t')