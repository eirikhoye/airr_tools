Github repository for the paper "T cell receptor repertoire sequencing reveals chemotherapy-driven clonal expansion in colorectal liver metastases" by Høye et al.

Adaptive ImmunoSEQ rearrangement files can be downloaded from:  
clients.adaptivebiotech.com  
email: hoye-review@adaptivebiotech.com  
password: hoye2022review

As well as from:

https://drive.google.com/drive/folders/1K0XESt0sMMNieD-YYv9GmFOOPs2PtHqB?usp=sharing


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


Overview of files required for the R_markdown analysis.

Data directory:
- SampleOverview_10-11-2021_1-47-39_PM.tsv: Sequencing parameters for all samples, including sequencing bathc 376276 and 12511620, downloaded from the immunoSEQ ANALYZER SampleOverview page. Loaded on line 57 in the markdown
- qcReport_kit_1.tsv: QC report from sequencing batch 376276, downloaded from immunoSEQ ANALYZER v2, samples with sequencing coverage < 5 where excluded from clonality analysis. Loaded on line 227 in the markdown
- qcReport_kit_2.tsv: QC report from sequencing batch 12511620, downloaded from immunoSEQ ANALYZER v2, samples with sequencing coverage < 5 where excluded from clonality analysis. Loaded on line 228 in the markdown
- metadata_all_v2.txt: metadata including COMET_ID and NACT group for all datasets. Loaded on line 60 of the markdown

hill_div directory contains tsv files output from the Hill_Diversity_v2.R script

envs directory contain yaml files for conda virtual environments needed to run scripts. Also needed for snakemake rules.

scripts directory contain the scripts neccessary for some precursor steps in the analysis. 
- Hill_Diversity_v2.R is used to make Hill diversity and evenness profiles. Script intended to be used with snakemake.
- Hill_Rarefy
