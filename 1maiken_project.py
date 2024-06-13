import argparse
import os
import pandas as pd
import numpy as np
import re

parser = argparse.ArgumentParser(description='Universal expression comparison\n', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metadata', required=False, help='TSV file with columns for sample ID, sex and tissue')
parser.add_argument('--matrix', required=False, help='File containing expression estimates')
args = parser.parse_args()

meta = args.metadata
matrix = args.matrix

ss= pd.read_csv(meta, sep="\t")
mat= pd.read_csv(matrix, sep="\t")
genecolumn= mat.columns[0]
mat = mat.set_index(genecolumn)

metadata_cats = ss.columns
print(metadata_cats)

ss_samplenames= list(ss["Sample"])
mat_samplenames= list(mat.columns)
assert len(np.unique(ss_samplenames)) == len(ss_samplenames)
assert len(np.unique(mat_samplenames)) == len(mat_samplenames)
assert set(ss_samplenames) == set(mat_samplenames)

def check_tissue(string):
    match = re.search(string, str(ss["Tissue"]), re.IGNORECASE)
    if match:
        return True
    else:
        return False

tissue_count=0
check_wholelfly= check_tissue("whole{0,2}.fly")
check_malpighian= check_tissue("malpighian{0,2}.tubule")
check_head= check_tissue("Head")
check_fatbody= check_tissue("fat{0,2}.body")
check_brain= check_tissue("brain")
for check in (check_malpighian, check_head, check_fatbody, check_brain):
    if check==True:
        tissue_count += 1
specificity = False
enrichment = False
if tissue_count >= 2:
    specificity = True
if tissue_count >=1 and check_wholelfly == 1:
    enrichment = True

check_female = False
check_male = False
for row in ss["Sex"]:
    if check_female == False:
        check_female= row.upper().startswith("F")
    if check_male == False:
        check_male = row.upper().startswith("M")
    if check_female == True and check_male == True:
        break

def check_age(age):
    match= re.search(age, str(ss['Age']), re.IGNORECASE)
    if match:
        return True
    else:
        return False
    
check_age_days = check_age(r'\d')
check_age_stage = check_age("[A-Z]")

print(str(ss["Age"]))