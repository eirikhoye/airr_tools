Github repository for the paper "T cell receptor repertoire sequencing reveals chemotherapy-driven clonal expansion in colorectal liver metastases" by Høye et al.

Adaptive ImmunoSEQ rearrangement files can be downloaded from:  
clients.adaptivebiotech.com  
email: hoye-review@adaptivebiotech.com  
password: hoye2022review

As well as from:

https://drive.google.com/drive/folders/1K0XESt0sMMNieD-YYv9GmFOOPs2PtHqB?usp=sharing


Overview of files required for the R_markdown analysis.
```
data/
- SampleOverview_10-11-2021_1-47-39_PM.tsv: Sequencing parameters for all samples, including sequencing bathc 376276 and 12511620, downloaded from the immunoSEQ ANALYZER SampleOverview page. Loaded on line 57 in the markdown
- qcReport_kit_1.tsv: QC report from sequencing batch 376276, downloaded from immunoSEQ ANALYZER v2, samples with sequencing coverage < 5 where excluded from clonality analysis. Loaded on line 227 in the markdown
- qcReport_kit_2.tsv: QC report from sequencing batch 12511620, downloaded from immunoSEQ ANALYZER v2, samples with sequencing coverage < 5 where excluded from clonality analysis. Loaded on line 228 in the markdown
- metadata_all_v2.txt: metadata including COMET_ID and NACT group for all datasets. Loaded on line 60 of the markdown

hill_div/ 
- directory contains tsv files output from the Hill_Diversity_v2.R script

envs/
- contain yaml files for conda virtual environments needed to run scripts. Also needed for snakemake rules.

scripts/
- contain the scripts neccessary for some precursor steps in the analysis. 
- Hill_Diversity.R is used to make Hill diversity and evenness profiles. Script part of snakemake rule.
- concat_Hill.R, cleans up diversity and evenness output into tidy dataframe.

```

To generate hill diversity profiles, edit the config.yaml file to contain the appropriate paths.

Because generating diversity profiles can be annoying, it is a good idea to run them through the snakemake workflow script, so each sample only has to be run once, regardless of whether additional samples are added to the analysis.

First, ensure that the rearrangement files are organized as one file per sample in the data/rearrangements directory. If you have all rearrangements merged into a single file, you can unmerge them in the same directory using: 
```
python scripts/Rearrangement_to_single_files.py </path/to/rearrangements_dir/> <rearrangements_filename.tsv>
```

To run the snakemake script, first install and activate a conda environment containing snakemake library, then use the following command:
```
# Install snakemake, follow:
# https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake

# Then, in the same directory as the SnakeMake file, run:
snakemake -j<number of cores> --configfile <path to config.yaml> --use-conda

# Note, it is important that the .yaml files are present in the envs/ directory, so that the neccessary conda environments are present.

# This will run the snakemake rules to generate diversity and evenness profiles in the results/<project_name> folder, and also concatenate them into a single tidy dataframe, as well as a dataframe with auc derived clonality values for each sample, which are easy to work with.

```

The hill_all.tsv can be loaded in R, and ggplot2 can easily produce publication quality figures with the following code:

```


```

![Screenshot]\(figures) 

















To set up imnet for network analysis, do:
```
# First create venv with 
virtualenv -p /usr/bin/python3.6 venv
source venv/bin/activate
pip install imnet pyspark
pip install findspark
pip install imnet

# must edit venv/lib/python3.6/site-packages/imnet/process_strings.py
# replace all instances of  idxs = range(nstrings)  to idxs = list(range(nstrings))

# First activate venv and make sure pyspark environment is set:
source venv/bin/activate
export PYSPARK_PYTHON=/home/jcdenton/imnet/venv/bin/python3.6
export PYSPARK_DRIVER_PYTHON=/home/jcdenton/imnet/venv/bin/python3.6

# run the test
pyhton test_imnet.py
```
