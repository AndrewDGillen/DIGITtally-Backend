import argparse
import pandas as pd
import re
from maiken_read_meta import *
import os
from sys import argv
import logging

def argparse_analysis():
    parser = argparse.ArgumentParser(description='Gene expression comparison\n', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--metadata', required=True, help='TSV file with columns for sample, sex, tissue and age')
    parser.add_argument('--matrix', required=True, help='File containing gene/sample matrix of gene expression estimates')
    parser.add_argument('-s', action='store_true', help='Run specificity True/False')
    parser.add_argument('--sthres', default = 1, type=float, help='Threshold ratio for specificity')
    parser.add_argument('-e', action='store_true', help='Run enrichment True/False')
    parser.add_argument('--ethres', default = 2.5, type=float, help='Threshold ratio for enrichment')
    parser.add_argument('--bthres', default = 2, type=float, help='Background threshold for gene count')
    parser.add_argument('--tissue', default = [], nargs='*', type = str, help='Tissue of interest')
    parser.add_argument('--permissive', default = [], nargs='*', type = str, help='Tissue to leave out of the analysis')
    parser.add_argument('--sex', default = [], nargs='*', type = str, help='Sex of interest')
    parser.add_argument('--age', default = [], nargs='*', type = str, help='Age of interest')
    parser.add_argument('--directory', required=True, type = str, help='Path to output directory')
    parser.add_argument('--logdir', required=True, type = str, help='Path to logging directory')
    parser.add_argument('--decimal', default='.', type = str, help='Character to recognize as decimal point in matrix values')
    
    return parser.parse_args()

class Tissue_expression_analyser:
    def __init__(self, arg_parser):
        self.directory = arg_parser.directory
        self.outputdir = self.check_outputdir()
        self.logdir = arg_parser.logdir 
        self.matrix = arg_parser.matrix
        self.meta = arg_parser.metadata
        self.decimal = arg_parser.decimal
        self.matrix_pd = self.read_matrix()
        self.meta_pd = read_meta(arg_parser)
        self.specificity = arg_parser.s
        self.enrichment = arg_parser.e
        self.sthres = float(arg_parser.sthres)
        self.ethres = float(arg_parser.ethres)
        self.bthres = float(arg_parser.bthres)
        self.tissue = arg_parser.tissue
        self.permissive = arg_parser.permissive
        self.sex = arg_parser.sex
        self.age = arg_parser.age
        
    def check_outputdir(self):
            output_dir= f"{self.directory}/Userdata analysis"
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
    
    def error_log(self, error):
        logger = logging.getLogger()
        if not logger.hasHandlers():
            logger.setLevel(logging.DEBUG)
            fh = logging.FileHandler(f'{self.logdir}/runtime_error.log')
            fh.setLevel(logging.ERROR)
            fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
            logger.addHandler(fh)
        return logger.error(error)
    
    def info_log(self, info):
        logger = logging.getLogger()
        if not logger.hasHandlers():
            logger.setLevel(logging.DEBUG)
            fh = logging.FileHandler(f'{self.logdir}/runtime_info.log')
            fh.setLevel(logging.INFO)
            fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
            logger.addHandler(fh)
        return logger.info(info)
    
    def read_matrix(self):
        try:
            assert self.matrix.endswith(".tsv")
            try: 
                matrix_pd = pd.read_csv(self.matrix, decimal=self.decimal, sep="\t", index_col=0, low_memory=False)
            except FileNotFoundError:
                self.error_log(f"Matrix file '{self.matrix}' does not exist. Analysis cancelled.")
                raise SystemExit(1)
        except AssertionError:
            self.error_log(f"Matrix file '{self.matrix}' needs to be in TSV format. Analysis cancelled.")
            raise SystemExit(1)
        matrix_pd = matrix_pd.apply(pd.to_numeric, errors="coerce")
        if matrix_pd.isnull().values.any():
            matrix_pd = matrix_pd.fillna(0)
            self.info_log(f"Matrix file '{self.matrix}' contains non-numeric datapoints, e.g. text, missing data, or invalid decimal point, which has been converted to 0.0")
        return matrix_pd

    def get_matrix_samples(self):
        samples= list(self.matrix_pd.columns)
        try:
            assert len(set(samples)) == len(samples)
            return samples
        except AssertionError:
            self.error_log(f"Matrix file '{self.matrix}' contains duplicate samples IDs.")
            raise SystemExit(1)
        
    def check_meta_matrix_samples(self):
        meta_samples = get_meta_samples(self.meta_pd)
        matrix_samples = self.get_matrix_samples()
        if set(meta_samples) == (set(matrix_samples)):
            return True
        elif set(meta_samples).issubset(set(matrix_samples)):
            self.info_log(f"Not all sample IDs in matrix '{self.matrix}' is present in metadata '{self.meta}'.")
            return True
        else:
            self.error_log(f"The sample IDs do not match between {self.meta} and {self.matrix}. Analysis cancelled.")
            raise SystemExit(1)

    def query_dict(self, tis):
        sex = self.sex
        age = self.age
        query = {}
        if bool(tis) == True:
            query['Tissue'] = [tis]
        if bool(sex) == True:
            query['Sex'] = sex
        if bool(age) == True:
            query['Age'] = age
        return query

    def filter_samples(self, query):
        filtered_metadata = self.meta_pd
        for key, values in query.items():
                filtered_metadata = filtered_metadata[filtered_metadata[key].isin(values)]
        return filtered_metadata['Sample']

    def run_specificity(self, tissue):
        if self.specificity:
            all_tissues = get_meta_tissues(self.meta_pd, False)
            other_targets= list(filter(lambda x: x!=tissue, self.tissue))
            non_target_means = {}
            for tis in all_tissues:
                if tis not in self.permissive and tis not in other_targets:
                    query = self.query_dict(tis)
                    query_samples = self.filter_samples(query)
                    if len(query_samples) < 3:
                        self.info_log(f"Specificity tissue '{tis}' only has {len(query_samples)} replicate(s) with age '{self.age}' and sex '{self.sex}'")
                    query_matrix = self.matrix_pd[query_samples]
                    query_mean = query_matrix.mean(axis=1)
                    query_mean= query_mean.mask(query_mean < self.bthres, self.bthres)
                    if tis == tissue:
                        target_mean = query_mean
                    else:
                        non_target_means[tis] = (list(query_mean))
            non_target_means = pd.DataFrame(non_target_means, index=self.matrix_pd.index).max(axis=1)
            result = target_mean/non_target_means
            result= result.apply(lambda x: True if x >= self.sthres else False)
            result.name = "Specificity"
            output = result[result == True].index
            with open(f"{self.outputdir}/{tissue}_specificity.tsv", "w") as specific_file:
                for gene in output:
                    specific_file.write(f"{gene}\n")
            return result

    def run_enrichment(self, tissue):
        if self.enrichment:
            all_tissues = get_meta_tissues(self.meta_pd, True)
            for tis in all_tissues:
                if find_wholefly(tis) == True or tis == tissue:
                    query = self.query_dict(tis)
                    query_samples = self.filter_samples(query)
                    if len(query_samples) < 3:
                        self.info_log(f"Enrichment tissue '{tis}' only has {len(query_samples)} replicate(s) with age '{self.age}' and sex '{self.sex}'")
                    query_matrix = self.matrix_pd[query_samples]
                    query_mean = query_matrix.mean(axis=1)
                    query_mean= query_mean.mask(query_mean < self.bthres, self.bthres)
                    if tis == tissue:
                        target_mean = query_mean
                    else:
                        whole_mean = query_mean
            result = target_mean/whole_mean
            result= result.apply(lambda x: True if x >= self.ethres else False)
            result.name = "Enrichment"
            output = result[result == True].index
            with open(f"{self.outputdir}/{tissue}_enrichment.tsv", "w") as enriched_file:
                for gene in output:
                    enriched_file.write(f"{gene}\n")
            return result
        
def execute_analysis(artificial_args=None):

    if artificial_args:
        arg_parser = Tissue_expression_analyser(artificial_args)
    else:
        arg_parser = Tissue_expression_analyser(argparse_analysis())

    if len(arg_parser.tissue) > 0:
        output_dir = arg_parser.outputdir
        arg_parser.check_meta_matrix_samples()
        combined_specificity = None
        combined_enrichment = None
        for tissue in arg_parser.tissue:
            specificity_result = arg_parser.run_specificity(tissue)
            enrichment_result = arg_parser.run_enrichment(tissue)
            if specificity_result is not None or enrichment_result is not None:
                merged = pd.concat([specificity_result, enrichment_result], axis=1)
                merged.to_csv(f"{output_dir}/{tissue}_expression.tsv", sep='\t')
            if combined_specificity is None:
                combined_specificity = specificity_result
            else:
                combined_specificity = combined_specificity & specificity_result
            if combined_enrichment is None:
                combined_enrichment = enrichment_result
            else:
                combined_enrichment = combined_enrichment & enrichment_result
        if specificity_result is not None or enrichment_result is not None:
            combined_merged = pd.concat([combined_specificity, combined_enrichment], axis=1)
            combined_merged.to_csv(f"{output_dir}/target_tissues_expression.tsv", sep='\t')
        if combined_specificity is not None:
            combined_specificity = combined_specificity[combined_specificity == True].index
            with open(f"{output_dir}/target_tissues_specificity.tsv", "w") as comb_spec_file:
                for gene in combined_specificity:
                    comb_spec_file.write(f"{gene}\n")
        if combined_enrichment is not None:
            combined_enrichment = combined_enrichment[combined_enrichment == True].index

            with open(f"{output_dir}/target_tissues_enrichment.tsv", "w") as comb_enrich_file:
                for gene in combined_enrichment:
                    comb_enrich_file.write(f"{gene}\n")
    arg_parser.info_log(f"Parameters used:\n\
    Directory: {arg_parser.directory}\n\
    Matrix file: {arg_parser.matrix}\n\
    Metadata file: {arg_parser.meta}\n\
    Decimal point: {arg_parser.decimal}\n\
    Run specificity: {arg_parser.specificity}\n\
    Run enrichment: {arg_parser.enrichment}\n\
    Specificity threshold:{arg_parser.sthres}\n\
    Enrichment threshold:{arg_parser.ethres}\n\
    Background threshold:{arg_parser.bthres}\n\
    Tissue(s) of interest: {", ".join(arg_parser.tissue)}\n\
    Permissive tissue(s): {", ".join(arg_parser.permissive)}\n\
    Sex(es) of interest: {", ".join(arg_parser.sex)}\n\
    Age(s) of interest: {", ".join(arg_parser.age)}")

    return output_dir

if __name__ == '__main__':
    execute_analysis()