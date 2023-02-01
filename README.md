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

Main directory:

SampleOverview_10-11-2021_1-47-39_PM.tsv: Sequencing parameters for all samples, including sequencing bathc 376276 and 12511620, downloaded from the immunoSEQ ANALYZER SampleOverview page. Loaded on line 57 in the markdown

rearrangement_details_kit_376276.tsv.zip: Processed sequencing reads from sequencing batch 376276, downloaded from the immunoSEQ ANALYZER Rearrangement_Details View page. Loaded on line 28 of the markdown

rearrangement_details_kit_12511620.tsv.zip: Processed sequencing reads from sequencing batch 1211620, downloaded from the immunoSEQ ANALYZER Rearrangement_Details View page. Loaded on line 34 of the markdown

qcReport_kit_1.tsv: QC report from sequencing batch 376276, downloaded from immunoSEQ ANALYZER v2, samples with sequencing coverage < 5 where excluded from clonality analysis. Loaded on line 227 in the markdown

qcReport_kit_2.tsv: QC report from sequencing batch 12511620, downloaded from immunoSEQ ANALYZER v2, samples with sequencing coverage < 5 where excluded from clonality analysis. Loaded on line 228 in the markdown

metadata_all_v2.txt: metadata including COMET_ID and NACT group for all datasets. Loaded on line 60 of the markdown

CombinedRearrangements_ALL_COMET.tsv: All sequencing reads for all datasets, used for morisita horn analysis, downloaded from immunoSEQ ANALYZER CombinedRearrangements page. Loaded on line 208 in the markdown




hill_div directory contains files that were output from the Hill_Diversity_v2.R script




miscellaneous folder with files that have been formatted slightly differently from the main metadata files because they are used by the R-markdown, but otherwise contain the same data as is present in metadata files in the main folder:

miscellaneous/SampleOverview_Kit-376276.tsv: Sequencing parameters for sequencing batch 376276, downloaded from the immunoSEQ ANALYZER SampleOverview View page. Loaded on line 46 of the R-markdown

miscellaneous/Overview_samples_kit_2_v2.tsv: file with information on input DNA concentration for sequencing batch 12511620. Loaded on line 48 of the R-markdown

miscellaneous/metadata_V2.txt: metadata including COMET_ID, NACT group for sequencing batch 376276, loaded on line 27 in the R-markdown.

miscellaneous/376276/metadata_v2.txt: loaded on line 221 in the R-markdown

miscellaneous/12511620//metadata_v2.txt: loaded on line 222 in the R-markdown
