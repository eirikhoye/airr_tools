import pandas as pd

files = snakemake.input

for file in files:
    pd.read_csv(file, sep='\t')
