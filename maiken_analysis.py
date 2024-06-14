import argparse
import pandas as pd
import re
from maiken_read_meta import *
import os

parser = argparse.ArgumentParser(description='Gene expression comparison\n', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-md', help='TSV file with columns for sample, sex, tissue and age')
parser.add_argument('-mx', help='File containing gene/sample matrix of gene expression estimates')
args = parser.parse_args()
matrix = args.mx
meta = args.md

def set_up_matrix(matrix_pd):
    genecolumn= matrix_pd.columns[0]
    matrix_pd = matrix_pd.set_index(genecolumn)
    return matrix_pd

def get_matrix_samples(matrix_pd):
    samples= list(matrix_pd.columns)
    return samples

def check_meta_matrix_samples(meta_pd, matrix_pd):
    meta_samples = get_meta_samples(meta_pd)
    matrix_samples = get_matrix_samples(matrix_pd)
    try:
        assert len(set(matrix_samples)) == len(matrix_samples)
        assert set(meta_samples) == set(matrix_samples)
    except AssertionError:
        print("The samples IDs either contain duplicates or do not match between the metadata and the matrix")

"""def make_database(meta_pd)
    if not os.path.isfile("metadata.db"):
        #A second file with the SQL create statements is opened to read.
        try:
            with open('SQL_DDL_2938235.txt', 'r') as sql_ddl:
                ddl_script = sql_ddl.read()
                db_manager.create_db(ddl_script)
        except FileNotFoundError:
            logger.error(f"SQL file 'SQL_DDL_2938235.txt' is not in the same directory as the program. Please ensure the file and program are in the working directory and try again.")
            raise SystemExit (1)
    else:
         logger.error(f'Database "{args.database}" already exists. Please try again with a different filename.')"""

def build_query_dict(tissue, sex, age):
    query = {}
    if bool(tissue) == True:
        query['Tissue'] == tissue
    if bool(sex) == True:
        query['Sex'] == sex
    if bool(age) == True:
        query['Age'] == age
    return query

def filter_samples(meta_pd, query):
    filtered_metadata = meta_pd
    for key, value in query.items():
        filtered_metadata = filtered_metadata[filtered_metadata[key] == value]
    return filtered_metadata['Sample']

def run_specificity(specificity, meta_pd, matrix_pd, tissue, sex, age):
    if specificity:
        
        query = build_query_dict(tissue, sex, age)
        query_samples = filter_samples(meta_pd, query)
        filtered_matrix = matrix_pd[query_samples]
        average_gene_expression = filtered_matrix.mean(axis=0)
        return average_gene_expression

meta= read_tsv(meta)
matrix= read_tsv(matrix)
meta_samples = get_meta_samples(meta)
matrix_samples = get_matrix_samples(matrix)
check_meta_matrix_samples(meta, matrix)
query_sex = "Female"
age_category = "Days"
query_age = 5
query_tissue = "Brain"
query_specificity = True
query_enrichment = True

#print(run_specificity.head(10))




