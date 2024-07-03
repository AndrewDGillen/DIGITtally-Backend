import argparse
import pandas as pd
import re
import logging

def meta_logger(log):
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(f'metadata_error.log')
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
    logger.addHandler(fh)
    return logger.error(log)

def argparse_meta():
    parser = argparse.ArgumentParser(description='Metadata reader\n', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--metadata', required=False, help='TSV file with columns for sample, sex, tissue and age')
    return parser.parse_args()

def read_meta(arg_parser):
    meta_file = arg_parser.metadata
    try:
        assert meta_file.endswith(".tsv")
        try:
            meta_pd = pd.read_csv(meta_file, sep="\t")
            meta_pd = meta_pd.astype(str)
            for column in meta_pd.columns:
                meta_pd[column]= meta_pd[column].str.strip()
                meta_pd[column]= meta_pd[column].str.replace(" ","_")
            return meta_pd
        except FileNotFoundError:
            meta_logger(f"Metadata file {meta_file} does not exist")
            raise SystemExit(1)
    except AssertionError:
        meta_logger(f"Metadata file needs to be in TSV format.")
        raise SystemExit(1)

def get_meta_cats(meta_pd):
    categories = set(meta_pd.columns)
    try:
        assert "Sample" in categories
        try:
            assert "Tissue" in categories
            return categories
        except AssertionError:
            meta_logger("Metadata does not contain tissue column")
            raise SystemExit(1)
    except AssertionError:
        meta_logger("Metadata does not contain sample column")
        raise SystemExit(1)

def get_meta_samples(meta_pd):
    samples= meta_pd["Sample"]
    if not len(set(samples)) == len(list(samples)):
        meta_logger("Metadata file contains duplicate samples. Analysis cancelled.")
        raise SystemExit(1)
    if samples.isnull().any():
        meta_logger("Metadata file is missing data in 'Sample' column")
        raise SystemExit(1)
    return list(samples)

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
            return set(no_whole)
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
    tissue_count= len(get_meta_tissues(meta_pd, False))
    if tissue_count > 1:
        return True
    else:
        return False

def check_enrichment(meta_pd):
    if check_wholelfly(meta_pd):
        tissue_count= len(get_meta_tissues(meta_pd, True))
        if tissue_count > 1:
            return True
    else:
        return False

def summary_dict():
    dict = {}
    meta= read_meta(argparse_meta())
    dict["meta_sexes"] = get_meta_sexes(meta)
    dict["age_category"] = get_age_categories(meta)
    dict["meta_ages"] = get_meta_ages(meta)
    dict["meta_tissues"] = get_meta_tissues(meta, False)
    dict["specificity"] = check_specificity(meta)
    dict["enrichment"] = check_enrichment(meta)
    return dict