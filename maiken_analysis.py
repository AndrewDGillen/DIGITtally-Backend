import argparse
import pandas as pd
import re
from maiken_read_meta import *
import os
from sys import argv

def argparse_analysis():
    parser = argparse.ArgumentParser(description='Gene expression comparison\n', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--metadata', required=True, help='TSV file with columns for sample, sex, tissue and age')
    parser.add_argument('--matrix', required=True, help='File containing gene/sample matrix of gene expression estimates')
    parser.add_argument('-s', action='store_true', help='Run specificity True/False')
    parser.add_argument('--sthres', default = 2, help='Threshold ratio for specificity')
    parser.add_argument('-e', action='store_true', help='Run enrichment True/False')
    parser.add_argument('--ethres', default = 2, help='Threshold ratio for enrichment')
    parser.add_argument('--bthres', default = 2, help='Background threshold for gene count')
    parser.add_argument('--tissue', required=True, help='Tissue of interest')
    parser.add_argument('--permissive', help='Tissue to leave out of the analysis')
    parser.add_argument('--sex', help='Sex of interest')
    parser.add_argument('--age', help='Age of interest')
    parser.add_argument('--directory', required=True, help='Path to output directory')
    parser.add_argument('--decimal', default='.', help='Character to recognize as decimal point')
    return parser.parse_args()


def read_matrix(arg_parser):
    matrix_file = arg_parser.matrix
    decimal=arg_parser.decimal
    try:
        assert matrix_file.endswith(".tsv")
        try: 
            matrix_pd = pd.read_csv(matrix_file, decimal=decimal, sep="\t", index_col=0, low_memory=False)
        except FileNotFoundError:
            print(f"Matrix file '{matrix_file}' does not exist")
            raise SystemExit(1)
    except AssertionError:
        print(f"Matrix file '{matrix_file}' needs to be in TSV format.")
        raise SystemExit(1)
    matrix_pd = matrix_pd.apply(pd.to_numeric, errors="coerce")
    if matrix_pd.isnull().values.any():
        matrix_pd = matrix_pd.fillna(0)
        print(f"Matrix file '{arg_parser.matrix}' contains non-numeric datapoints, e.g. text or invalid decimal point, which has been converted to 0.0")
    #matrix_pd.astype(float)
    """try:
        matrix_pd = matrix_pd.apply(pd.to_numeric).fillna(0)
    except ValueError:
        print(f"Matrix file '{arg_parser.matrix}' contains non-numeric values")
        raise SystemExit(1)"""
    matrix_samples = get_matrix_samples(matrix_pd)
    try:
        assert len(set(matrix_samples)) == len(matrix_samples)
    except AssertionError:
        print(f"Matrix file '{arg_parser.matrix}' contains duplicate samples IDs")
        raise SystemExit(1)
    return matrix_pd

def get_matrix_samples(matrix_pd):
    samples= list(matrix_pd.columns)
    return samples

def check_meta_matrix_samples(arg_parser):
    meta_pd = read_meta(arg_parser)
    matrix_pd = read_matrix(arg_parser)
    meta_samples = get_meta_samples(meta_pd)
    matrix_samples = get_matrix_samples(matrix_pd)
    if set(meta_samples) == (set(matrix_samples)):
        return True
    elif set(meta_samples).issubset(set(matrix_samples)):
        print(f"INFO: Not all sample IDs in matrix '{arg_parser.matrix}' is present in metadata '{arg_parser.metadata}'. Analysis continued.")
        return True
    else:
        print(f"The samples IDs do not match between {arg_parser.metadata} and {arg_parser.matrix}. Analysis cancelled.")
        raise SystemExit(1)

def query_dict(arg_parser, tis):
    sex = arg_parser.sex
    age = arg_parser.age
    query = {}
    if bool(tis) == True:
        query['Tissue'] = tis
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

def run_specificity(arg_parser, output_dir):
    specificity = arg_parser.s
    if specificity:
        meta_pd = read_meta(arg_parser)
        matrix_pd = read_matrix(arg_parser)
        threshold = float(arg_parser.sthres)
        tissue = arg_parser.tissue
        permissive = arg_parser.permissive
        all_tissues = get_meta_tissues(meta_pd, False)
        non_target_means = {}
        for tis in all_tissues:
            if tis is not permissive:
                query = query_dict(arg_parser, tis)
                query_samples = filter_samples(meta_pd, query)
                if len(query_samples) < 3:
                    print(f"Specificity tissue '{tis}' only has {len(query_samples)} replicate(s) with age '{arg_parser.age}' and sex '{arg_parser.sex}'")
                query_matrix = matrix_pd[query_samples]
                query_mean = query_matrix.mean(axis=1)
                if tis == tissue:
                    main_mean = query_mean
                else:
                    non_target_means[tis] = (list(query_mean))
        non_target_means = pd.DataFrame(non_target_means, index=matrix_pd.index).max(axis=1)
        result = main_mean/non_target_means
        result= result.apply(lambda x: True if x >= threshold else False)
        result.name = "Specificity"
        output = result[result == True].index
        with open(f"{output_dir}/Specific_to_{tissue}.tsv", "w") as specific_file:
            for gene in output:
                specific_file.write(f"{gene}\n")
        return result

def run_enrichment(arg_parser, output_dir):
    enrichment = arg_parser.e
    if enrichment:
        meta_pd = read_meta(arg_parser)
        matrix_pd = read_matrix(arg_parser)
        threshold = float(arg_parser.ethres)
        background_threshold = float(arg_parser.bthres)
        tissue = arg_parser.tissue
        all_tissues = get_meta_tissues(meta_pd, True)
        for tis in all_tissues:
            if find_wholefly(tis) == True or tis == tissue:
                query = query_dict(arg_parser, tis)
                query_samples = filter_samples(meta_pd, query)
                if len(query_samples) < 3:
                    print(f"Enrichment tissue '{tis}' only has {len(query_samples)} replicate(s) with age '{arg_parser.age}' and sex '{arg_parser.sex}'")
                query_matrix = matrix_pd[query_samples]
                query_mean = query_matrix.mean(axis=1)
                query_mean= query_mean.mask(query_mean < background_threshold, 2)
                if tis == tissue:
                    main_mean = query_mean
                else:
                    whole_mean = query_mean
        result = main_mean/whole_mean
        result= result.apply(lambda x: True if x >= threshold else False)
        result.name = "Enrichment"
        output = result[result == True].index
        with open(f"{output_dir}/Enriched_in_{tissue}.tsv", "w") as enriched_file:
            for gene in output:
                enriched_file.write(f"{gene}\n")
        return result

def check_outputdir(output_dir):
    output_dir= f"{output_dir}/Expression analysis"
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)
    else:
        extension = 1
        intermediate = output_dir
        while os.path.exists(output_dir):
            extension += 1
            output_dir = f"{intermediate}{extension}"
        os.mkdir(output_dir)
    return output_dir
        
def execute_analysis():
    arg_parser = argparse_analysis()
    tissue = arg_parser.tissue
    check_meta_matrix_samples(arg_parser)
    output_dir = check_outputdir(arg_parser.directory)
    specificity_result = run_specificity(arg_parser, output_dir)
    enrichment_result = run_enrichment(arg_parser, output_dir)
    merged = pd.concat([specificity_result, enrichment_result], axis=1)
    merged.to_csv(f"{output_dir}/Genes_in_{tissue}.tsv", sep='\t')

    with open(F"{output_dir}/Query_parameters.tsv", "w") as param_file:
        param_file.write(f"Commandline used:\n")
        for param in argv[1:]:
            param_file.write(f"{param} ")



