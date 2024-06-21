import argparse
import pandas as pd
import re
from maiken_read_meta import *
import os

def argparse_analysis():
    parser = argparse.ArgumentParser(description='Gene expression comparison\n', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--metadata', help='TSV file with columns for sample, sex, tissue and age')
    parser.add_argument('--matrix', help='File containing gene/sample matrix of gene expression estimates')
    parser.add_argument('-s', action='store_true', help='Run specificity True/False')
    parser.add_argument('--sthres', default = 2, help='Threshold ratio for specificity')
    parser.add_argument('-e', action='store_true', help='Run enrichment True/False')
    parser.add_argument('--ethres', default = 2, help='Threshold ratio for enrichment')
    parser.add_argument('--tissue', required=True, help='Tissue of interest')
    parser.add_argument('--sex', help='Sex of interest')
    parser.add_argument('--age', help='Age of interest')
    parser.add_argument('--directory', help='Path to output directory')
    #return parser


def read_matrix(arg_parser):
    parser = arg_parser()
    matrix_pd = parser.matrix
    matrix_pd = pd.read_csv(matrix_pd, sep="\t")
    genecolumn = matrix_pd.columns[0]
    matrix_pd = matrix_pd.set_index(genecolumn).apply(pd.to_numeric)
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

def query_dict(arg_parser):
    parser = arg_parser()
    tissue = parser.tissue
    sex = parser.sex
    age = parser.age
    query = {}
    if bool(tissue) == True:
        query['Tissue'] = tissue
    if bool(sex) == True:
        query['Sex'] = sex
    if bool(age) == True:
        query['Age'] = age
    return query

def filter_samples(meta_pd, query):
    filtered_metadata = meta_pd
    for key, value in query.items():
        filtered_metadata = filtered_metadata[filtered_metadata[key] == value]
    return filtered_metadata['Sample']

def write_file(data, filename, output_folder):
    """output_folder= f"{output_folder}/Expression analysis"
    if os.path.exists(output_folder) == False:
        os.mkdir(output_folder)
    else:
        extension = 1
        while os.path.exists(output_folder):
            extension += 1
            output_folder = f"{output_folder}{extension}"
        os.mkdir(output_folder)"""
    data.to_csv(f"{output_folder}/{filename}.tsv", sep='\t')

def run_specificity(arg_parser):
    parser = arg_parser()
    specificity = parser.s
    meta_pd = parser.metadata
    matrix_pd = parser.matrix
    threshold = parser.sthres
    tissue = parser.tissue
    output_folder = parser.directory
    if specificity:
        all_tissues = get_meta_tissues(meta, False)
        non_target_means = {}
        tissue = query_dict(arg_parser)
        tissue = tissue['Tissue']
        for tis in all_tissues:
            query = query_dict(arg_parser)
            query_samples = filter_samples(meta_pd, query)
            query_matrix = matrix_pd[query_samples]
            query_mean_expression = query_matrix.mean(axis=1)
            if tis == tissue:
                main_mean = query_mean_expression
            else:
                non_target_means[tis] = (list(query_mean_expression))
        non_target_means = pd.DataFrame(non_target_means, index=matrix_pd.index).max(axis=1)
        result = main_mean/non_target_means
        result= result.apply(lambda x: True if x >= threshold else False)
        result.name = "Specificity"
        output = result[result == True]
        write_file(output, f"Specific_to_{tissue}", output_folder)
        return result

def run_enrichment(arg_parser):
    parser = arg_parser()
    enrichment = parser.e
    meta_pd = parser.metadata
    matrix_pd = parser.matrix
    threshold = parser.ethres
    tissue = parser.tissue
    output_folder = parser.directory
    if enrichment:
        all_tissues = get_meta_tissues(meta, True)
        tissue = query_dict(parser)
        tissue = tissue['Tissue']
        for tis in all_tissues:
            if find_wholefly(tis) == True or tis == tissue:
                query = query_dict(parser)
                query_samples = filter_samples(meta_pd, query)
                query_matrix = matrix_pd[query_samples]
                query_mean_expression = query_matrix.mean(axis=1)
                if tis == tissue:
                    main_mean = query_mean_expression
                else:
                    whole_mean = query_mean_expression
        result = main_mean/whole_mean
        result= result.apply(lambda x: True if x >= threshold else False)
        result.name = "Enrichment"
        output = result[result == True]
        write_file(output, f"Enriched_in_{tissue}", output_folder)
        return result

def check_outputdir(output_folder):
    output_folder= f"{output_folder}/Expression analysis"
    if os.path.exists(output_folder) == False:
        os.mkdir(output_folder)
    else:
        extension = 1
        intermediate = output_folder
        while os.path.exists(output_folder):
            extension += 1
            output_folder = f"{intermediate}{extension}"
        os.mkdir(output_folder)
    return output_folder
        
def final_output(arg_parser):
    parser = arg_parser()
    tissue = parser.tissue
    output_folder = parser.directory
    output_folder = check_outputdir(output_folder)
    specificity_result = run_specificity(parser)
    enrichment_result = run_enrichment(parser)
    merged = pd.concat([specificity_result, enrichment_result], axis=1)
    write_file(merged, f"All_genes_{tissue}.tsv", output_folder)
    #Make query file

"""specificity_threshold = 2
enrichment_threshold = 2
meta_samples = get_meta_samples(meta)
matrix_samples = get_matrix_samples(matrix)
query_sex = False
age_category = "Days"
query_age = False
query_tissue = "Head"
query_specificity = True
query_enrichment = True
output_path= "C:/Users/maike/Documents/University/Masters/Project"""
meta = read_meta(argparse_analysis())
matrix = read_matrix(argparse_analysis())
check_meta_matrix_samples(meta, matrix)
final_output(argparse_analysis)




