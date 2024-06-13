import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description='Gene expression comparison\n', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metadata', required=False, help='TSV file with columns for sample, sex, tissue and age')
args = parser.parse_args()
meta = args.metadata

def read_tsv(tsv):
    frame = pd.read_csv(tsv, sep="\t")
    return frame

def get_meta_cats(meta_pd):
    categories = set(meta_pd.columns)
    assert "Sample" in categories
    return categories

def get_meta_samples(meta_pd):
    samples= list(meta_pd["Sample"])
    assert len(set(samples)) == len(samples)
    return samples

def get_age_categories(meta_pd):
    if "Age" in get_meta_cats(meta_pd):  
        def search_age(pattern, meta_pd):
            match= re.search(pattern, str(set(meta_pd['Age'])), re.IGNORECASE)
            if match:
                return True
            else:
                return False
        days =search_age(r'\d', meta_pd)
        stage =search_age("[A-Za-z]", meta_pd)
        if days == True and stage == True:
            return("Mixed")
        elif days == True:
            return("Days")
        elif stage == True:
            return("Stage")
    else:
        return None
    
def get_meta_ages(meta_pd):
    if "Age" in get_meta_cats(meta_pd):  
        ages = set(meta_pd["Age"])
        return ages
    else:
        return None   

def get_meta_sexes(meta_pd):
    if "Sex" in get_meta_cats(meta_pd):
        sexes = set(meta_pd['Sex'])
        return sexes
    else:
        return None
    
def get_meta_tissues(meta_pd):
    if "Tissue" in get_meta_cats(meta_pd):
        tissues = set(meta_pd['Tissue'])
        return tissues
    else:
        return None

def check_wholelfly(meta_pd):
    if "Tissue" in get_meta_cats(meta_pd):
        tissues= get_meta_tissues(meta_pd)
        wholefly = bool(re.search("whole", str(tissues), re.IGNORECASE))
        if wholefly:
            return True
        else:
            return False
    else:
        return False

def check_specificity(meta_pd):
    tissues= get_meta_tissues(meta_pd)
    tissue_count=len(tissues)
    if check_wholelfly(meta_pd):
        if tissue_count > 2:
            return True
    elif tissue_count > 1:
        return True
    else:
        return False

def check_enrichment(meta_pd):
    if check_wholelfly(meta_pd):
        tissues= get_meta_tissues(meta_pd)
        tissue_count=len(tissues)
        if tissue_count > 1:
            return True
    else:
        return False





