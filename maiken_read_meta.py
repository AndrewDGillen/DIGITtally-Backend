import argparse
import pandas as pd
import re

def argparse_meta():
    parser = argparse.ArgumentParser(description='Metadata reader\n', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--metadata', required=False, help='TSV file with columns for sample, sex, tissue and age')
    return parser.parse_args()

def read_meta(arg_parser):
    meta_pd = arg_parser.metadata
    meta_pd = pd.read_csv(meta_pd, sep="\t")
    return meta_pd

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
    
def find_wholefly(string):
    whole = bool(re.search("whole", string, re.IGNORECASE))
    return whole    

def get_meta_tissues(meta_pd, whole):
    if "Tissue" in get_meta_cats(meta_pd):
        tissues = set(meta_pd['Tissue'])
        if whole == False:
            no_whole = []
            for tis in tissues:
                if find_wholefly(tis) == False:
                    no_whole.append(tis)
            return no_whole
        else:
            return tissues
    else:
        return None

def check_wholelfly(meta_pd):
    if "Tissue" in get_meta_cats(meta_pd):
        tissues= get_meta_tissues(meta_pd, True)
        wholefly = find_wholefly(str(tissues))
        return wholefly
    else:
        return False

def check_specificity(meta_pd):
    tissues= get_meta_tissues(meta_pd, False)
    tissue_count=len(tissues)
    if tissue_count > 1:
        return True
    else:
        return False

def check_enrichment(meta_pd):
    if check_wholelfly(meta_pd):
        tissues= get_meta_tissues(meta_pd, True)
        tissue_count=len(tissues)
        if tissue_count > 1:
            return True
    else:
        return False

def summary_dict():
    dict = {}
    meta= read_meta(argparse_meta)
    dict["meta_sexes"] = get_meta_sexes(meta)
    dict["age_category"] = get_age_categories(meta)
    dict["meta_ages"] = get_meta_ages(meta)
    dict["meta_tissues"] = get_meta_tissues(meta, False)
    dict["specificity"] = check_specificity(meta)
    dict["enrichment"] = check_enrichment(meta)
    return dict