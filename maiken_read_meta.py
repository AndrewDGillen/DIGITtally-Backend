import argparse
import pandas as pd
import re
import logging
import json
import os 

config_file = '/media/ag_jdowlab/Seagate/DIGITtally-Staging/dt-config.json'
if os.path.isfile(config_file):
    with open(config_file) as config_file_in:
        config = json.load(config_file_in)
else:
    config={'LogLocation':os.getcwd()}

def meta_logger(log):
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(f'{config["LogLocation"]}/metadata_error.log')
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
    logger.addHandler(fh)
    return logger.error(log)

def argparse_meta():
    parser = argparse.ArgumentParser(description='Metadata reader\n', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--metadata', help='TSV file with columns for sample, sex, tissue and age')
    return parser.parse_args()

def read_meta(arg_parser):
    """Parses the metadata tsv file into a pandas dataframe"""
    meta_file = arg_parser.metadata
    try:
        assert meta_file.endswith(".tsv")
        try:
            meta_pd = pd.read_csv(meta_file, sep="\t")
            #The pandas dataframe is returned
            return meta_pd
        except FileNotFoundError:
            meta_logger(f"Metadata file {meta_file} does not exist")
            raise SystemExit(1)
    except AssertionError:
        meta_logger(f"Metadata file needs to be in TSV format.")
        raise SystemExit(1)

#File format and presence of required columns is checked
def get_meta_cats(meta_pd):
    """Extract the column headers (metadata categories) from the metadata dataframe"""
    categories = list(meta_pd.columns)
    try:
        assert len(categories) == len(set(categories))
        try:
            assert "Sample" in categories
            try:
                assert "Tissue" in categories
                return set(categories)
            except AssertionError:
                meta_logger("Metadata does not contain 'Tissue' column")
                raise SystemExit(1)
        except AssertionError:
            meta_logger("Metadata does not contain 'Sample' column")
            raise SystemExit(1)
    except AssertionError:
        meta_logger("Metadata contains duplicate columns")
        raise SystemExit(1)

#Checking for duplicates and empty datapoints 
def get_meta_samples(meta_pd):
    """Extract the sample IDs from the metadata"""
    samples= meta_pd["Sample"]
    if not len(set(samples)) == len(list(samples)):
        meta_logger("Metadata file contains duplicate sample IDs.")
        raise SystemExit(1)
    if samples.isnull().any():
        meta_logger("Metadata file is missing data in 'Sample' column")
        raise SystemExit(1)
    #Keeping the sample IDs stringtype ensures the samples of metadata and matrix can be compared, in case numerical IDs are rendered differently in each script
    samples = samples.astype(str)
    return list(samples)

#This is for user information and design of the DIGITtally options pages
def get_age_categories(meta_pd):
    """Identify whether the ages are given in days, life stages, or a mix"""
    if "Age" in get_meta_cats(meta_pd):  
        def search_age(pattern):
            """Regex search in the 'Age' column of metadata"""
            #The entire age column is turned into a string that can be searched for numbers and letters
            match= re.search(pattern, str(set(meta_pd['Age'].dropna())), re.IGNORECASE)
            if match:
                return True
            else:
                return False
        #Numbers and letters are searched for separately
        days = search_age(r'\d')
        stage = search_age("[A-Za-z]")
        if days == True and stage == True:
            return("Mixed")
        elif days == True:
            return("Days")
        elif stage == True:
            return("Stage")
    else:
        return None

def get_meta_ages(meta_pd):
    """Extract all the ages given in the metadata"""
    if "Age" in get_meta_cats(meta_pd):
        #The ages are kept as stringtype as they are only treated categorically and never as numericals
        ages = set(meta_pd["Age"].dropna().astype(str))
        return ages
    else:
        return None   

def get_meta_sexes(meta_pd):
    """Extract all the sexes given in the metadata"""
    if "Sex" in get_meta_cats(meta_pd):
        #The sexes are kept as stringtype in case the user uses numericals in place of sexes
        sexes = set(meta_pd['Sex'].dropna().astype(str))
        return sexes
    else:
        return None

#This function will be used to identify if a tissue should be treated as "wholefly" tissues
def find_wholefly(string):
    '''Check if a string contains the word "whole" (case insensitive)'''
    whole = bool(re.search("whole", string, re.IGNORECASE))
    return whole    

def get_meta_tissues(meta_pd, whole):
    '''Extract a set of the tissues in the metadata'''
    if "Tissue" in get_meta_cats(meta_pd):
        tissues = set(meta_pd['Tissue'].dropna().astype(str))
        whole_tissues = 0
        no_whole = []
        for tis in tissues:
            #Checking if there is more than one wholefly "tissue", as the enrichment analysis requires only one wholefly "tissue"
            if find_wholefly(tis) == True:
                whole_tissues += 1
                if whole_tissues > 1:
                    meta_logger('More than one tissue contain the word "whole" (case insensitive)')
                    raise SystemExit(1)
            #A list is made without the wholefly tissue, for the user to select target tissues from, and to loop over in the specificity analysis
            else:
                no_whole.append(tis.strip())
        #Sets of tissues with or without the wholefly "tissue" can be returned depending on need
        if whole == False:
            return set(no_whole)
        else:
            return set(x.strip() for x in tissues)
    else:
        return None

def check_wholelfly(meta_pd):
    '''Check presence of any wholefly "tissue" in the metadata'''
    if "Tissue" in get_meta_cats(meta_pd):
        tissues= get_meta_tissues(meta_pd, True)
        wholefly = find_wholefly(str(tissues))
        return wholefly
    else:
        return False

#The specificity analysis is only present if there's at least 2 non-wholefly tissues present
def check_specificity(meta_pd):
    '''Check whether a specificity analysis is possible based on the tissues present'''
    #Finding the number of non-wholefly tissues
    tissue_count= len(get_meta_tissues(meta_pd, False))
    if tissue_count > 1:
        return True
    else:
        return False

#The enrichment analysis is only possible if there's a wholefly "tissue" and at least 1 non-wholefly tissue present
def check_enrichment(meta_pd):
    '''Check whether an enrichment analysis is possible based on the tissues present'''
    if check_wholelfly(meta_pd):
        #Finding the number of tissues present including the wholefly "tissue"
        tissue_count= len(get_meta_tissues(meta_pd, True))
        if tissue_count > 1:
            return True
    else:
        return False

#This function acts as a main function that runs all the others, as well as returning a summary of the metadata content as a dictionary
def summary_dict():
    dict = {}
    meta= read_meta(argparse_meta())
    #The samples function needs to runs purely for error handling
    get_meta_samples(meta)
    dict["meta_sexes"] = get_meta_sexes(meta)
    dict["age_category"] = get_age_categories(meta)
    dict["meta_ages"] = get_meta_ages(meta)
    dict["meta_tissues"] = get_meta_tissues(meta, False)
    dict["specificity"] = check_specificity(meta)
    dict["enrichment"] = check_enrichment(meta)
    return dict

#This an alternative function to run the script which takes an argparse.Namespace object as an argumant, which parses the args from another script as if it was the expected commandline argparser
def get_meta_settings(metafile):
    meta_dict = {}
    meta= read_meta(metafile)
    get_meta_samples(meta)
    meta_dict["meta_sexes"] = get_meta_sexes(meta)
    meta_dict["age_category"] = get_age_categories(meta)
    meta_dict["meta_ages"] = get_meta_ages(meta)
    meta_dict["meta_tissues"] = get_meta_tissues(meta, False)
    meta_dict["specificity"] = check_specificity(meta)
    meta_dict["enrichment"] = check_enrichment(meta)
    return meta_dict
