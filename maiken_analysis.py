import argparse
import pandas as pd
import re
from maiken_read_meta import *
import os
from sys import argv
import logging

#Setting everything, including the argparser, up as a function helps organising the code and to import sections into the main script.
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

#Making a class to initiate and amend (e.g. object type and conversion into pandas dataframes) the args just once
#Including the amended args in the __init__ prevents errors throughout the script, such as forgetting to set the type at every use
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

    #This function creates a directory at users chosen destination but alters the directory name in case it already exists
    def check_outputdir(self):
        """Create uniquely names directory for output files""" 
        output_dir= f"{self.directory}/Userdata analysis"
        if os.path.exists(output_dir) == False:
            os.mkdir(output_dir)
        else:
            extension = 2
            #The baseline directory is stored as an intermediate to allow the extension to count up rather than appending a digit
            intermediate = output_dir
            while os.path.exists(output_dir):
                extension += 1
                output_dir = f"{intermediate}{extension}"
            os.mkdir(output_dir)
        return output_dir
    
    #Separate logs for errors and info are set up to allow easy lookup of reason for termination (all errors result in systemexit)
    def error_log(self, error):
        """Log error level messages"""
        logger = logging.getLogger()
        if not logger.hasHandlers():
            logger.setLevel(logging.DEBUG)
            fh = logging.FileHandler(f'{self.logdir}/runtime_error.log')
            fh.setLevel(logging.ERROR)
            fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
            logger.addHandler(fh)
        return logger.error(error)
    
    #The info log is intended to contain messages about the uploaded data, e.g. in case of missing values, low sample numbers, etc.
    #No errors preventing the script to run will be logged in the info log
    def info_log(self, info):
        """Log info level messages"""
        logger = logging.getLogger()
        if not logger.hasHandlers():
            logger.setLevel(logging.DEBUG)
            fh = logging.FileHandler(f'{self.logdir}/runtime_info.log')
            fh.setLevel(logging.INFO)
            fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
            logger.addHandler(fh)
        return logger.info(info)
    
    #This function reads in the tsv as a pandas dataframe and sets the indexes and types as needed for further analysis
    def read_matrix(self):
        """Read in a csv as a pandas dataframe with float type data"""
        try:
            assert self.matrix.endswith(".tsv")
            try: 
                #low_memory is set to False in case of any strings appearing in the dataframe
                matrix_pd = pd.read_csv(self.matrix, decimal=self.decimal, sep="\t", index_col=0, low_memory=False)
            except FileNotFoundError:
                self.error_log(f"Matrix file '{self.matrix}' does not exist.")
                raise SystemExit(1)
        except AssertionError:
            self.error_log(f"Matrix file '{self.matrix}' needs to be in TSV format.")
            raise SystemExit(1)
        #Using the coerce argument converts any non-numercial data into NA floats, which works for info logging and can be converted into value 0 (float) along with other empty data input
        matrix_pd = matrix_pd.apply(pd.to_numeric, errors="coerce")
        if matrix_pd.isnull().values.any():
            matrix_pd = matrix_pd.fillna(0)
            self.info_log(f"Matrix file '{self.matrix}' contains non-numeric datapoints, e.g. text, missing data, or invalid decimal point, which has been converted to 0.0")
        return matrix_pd

    #Checking for no duplicate sample IDs prevents errors when extracting specific columns during the analysis
    def get_matrix_samples(self):
        """Extract sample IDs from the column headers of the matrix and check there are no duplicates"""
        samples= list(self.matrix_pd.columns)
        try:
            assert len(set(samples)) == len(samples)
            return samples
        except AssertionError:
            self.error_log(f"Matrix file '{self.matrix}' contains duplicate samples IDs.")
            raise SystemExit(1)
    
    #Checking that the sample IDs match between the two files also prevents errors when extracting columns from the matrix during the analysis
    def check_meta_matrix_samples(self):
        """Compare sample IDs present in metadata and matrix"""
        meta_samples = get_meta_samples(self.meta_pd)
        matrix_samples = self.get_matrix_samples()
        if set(meta_samples) == (set(matrix_samples)):
            return True
        #Continueing the analysis if the metadata sample IDs are a subset of those of the matrix allows the user to work with different subsets of samples in the matrix at a time
        #This might be useful for example if a user has different experimental conditions where only using the tissue, age and sex parameters to filter the samples is not sufficient
        elif set(meta_samples).issubset(set(matrix_samples)):
            self.info_log(f"Not all sample IDs in matrix '{self.matrix}' is present in metadata '{self.meta}'.")
            return True
        else:
            self.error_log(f"The sample IDs do not match between {self.meta} and {self.matrix}.")
            raise SystemExit(1)

    #Using a dictionary with keys of the same name as the metadata columns makes it easy to loop over the combinations in case of multiple ages/sexes of interest
    def query_dict(self, tis):
        """Create dictionary to filter the metadata by"""
        sex = self.sex
        age = self.age
        query = {}
        if bool(tis) == True:
            query['Tissue'] = [tis]
        if bool(sex) == True:
            query['Sex'] = self.sex
        if bool(age) == True:
            query['Age'] = self.age
        return query

    def filter_samples(self, query):
        """Filter the metadata by using a dictionary"""
        filtered_metadata = self.meta_pd
        for key, values in query.items():
                filtered_metadata = filtered_metadata[filtered_metadata[key].isin(values)]
        return filtered_metadata['Sample']

    #This function runs the specificity analysis for a single target tissue (specified as an argument) using the query and filtering functions
    def run_specificity(self, tissue):
        """Run the specificity analysis where target-tissue expression is compared to non-target-tissue expression"""
        if self.specificity:
            all_tissues = get_meta_tissues(self.meta_pd, False)
            #By the specificity definition, a gene can only be specific to one tissue.
            #A list object with the remaining target tissues is made to leave these out of the comparison
            #In case of e.g. 3 target tissues, these tissues simply need to be the top 3 tissues for a gene, rather than all 3 each having the highest expression (impossible)
            other_targets= list(filter(lambda x: x!=tissue, self.tissue))
            #A dictionary is made to store non-target means while looping, as a dictionary easily can be converted into a dataframe
            non_target_means = {}
            for tis in all_tissues:
                if tis not in self.permissive and tis not in other_targets:
                    query = self.query_dict(tis)
                    query_samples = self.filter_samples(query)
                    #A log is made in case of low sample size
                    if len(query_samples) < 3:
                        self.info_log(f"Specificity tissue '{tis}' only has {len(query_samples)} replicate(s) with age '{self.age}' and sex '{self.sex}'")
                    #The samples fulfilling the metadata query are subsetted in the matrix
                    query_matrix = self.matrix_pd[query_samples]
                    query_mean = query_matrix.mean(axis=1)
                    #Any mean lower than the accepted RNAseq detection baseline is increased to the defined baseline
                    query_mean= query_mean.mask(query_mean < self.bthres, self.bthres)
                    if tis == tissue:
                        target_mean = query_mean
                    else:
                        non_target_means[tis] = (list(query_mean))
            #Only the highest mean expression per gene is kept, as only this is relevant to compare the target mean to
            non_target_means = pd.DataFrame(non_target_means, index=self.matrix_pd.index).max(axis=1)
            result = target_mean/non_target_means
            result= result.apply(lambda x: True if x >= self.sthres else False)
            result.name = "Specificity"
            #A file is created with a list of all genes found to be specific to the current target tissue
            output = result[result == True].index
            with open(f"{self.outputdir}/{tissue}_specificity.tsv", "w") as specific_file:
                for gene in output:
                    specific_file.write(f"{gene}\n")
            #The boolean values for all genes is returned at the end of the function for writing more files and finding genes specific to all target tissues.
            return result

    #This function runs the enrichment analysis for a single target tissue (specified as an argument) using the query and filtering functions
    def run_enrichment(self, tissue):
        """Run the enrichment analysis where target-tissue expression is compared to whole organism expression"""
        if self.enrichment:
            all_tissues = get_meta_tissues(self.meta_pd, True)
            for tis in all_tissues:
                #The mean expression will only be calculated for the target tissue and the wholefly "tissue"
                if find_wholefly(tis) == True or tis == tissue:
                    query = self.query_dict(tis)
                    query_samples = self.filter_samples(query)
                    #A log is made in case of low sample size
                    if len(query_samples) < 3:
                        self.info_log(f"Enrichment tissue '{tis}' only has {len(query_samples)} replicate(s) with age '{self.age}' and sex '{self.sex}'")
                    #The samples fulfilling the metadata query are subsetted in the matrix
                    query_matrix = self.matrix_pd[query_samples]
                    query_mean = query_matrix.mean(axis=1)
                    #Any mean lower than the accepted RNAseq detection baseline is increased to the defined baseline
                    query_mean= query_mean.mask(query_mean < self.bthres, self.bthres)
                    if tis == tissue:
                        target_mean = query_mean
                    else:
                        #The metadata script contains error handling to ensure only one tissue present is treated as the wholefly tissue
                        #That prevents the whole_mean object from being overwritten 
                        whole_mean = query_mean
            result = target_mean/whole_mean
            result= result.apply(lambda x: True if x >= self.ethres else False)
            result.name = "Enrichment"
            #A file is created with a list of all genes found to be enriched in the current target tissue
            output = result[result == True].index
            with open(f"{self.outputdir}/{tissue}_enrichment.tsv", "w") as enriched_file:
                for gene in output:
                    enriched_file.write(f"{gene}\n")
            #The boolean values for all genes is returned at the end of the function for writing more files and finding genes enriched in all target tissues.
            return result

#This function combines all the previous functions into one big analysis with background checks. 
#This makes it easy to run the script from another main script (i.e. the DIGITtally backbone)
def execute_analysis(artificial_args=None):
    """Execute the specificity and enrichment analysis including output directory checks and compatibility between the matrix the metadata"""
    
    #This if statement (combined with the "artificial_args=None" argument) allows the script to run either from commanline using the "argparse_analysis" or
    #Run the funtion using an argparse.Namespace object as an argumant, which parses the args from another script as if it was a commandline argparser
    if artificial_args:
        arg_parser = Tissue_expression_analyser(artificial_args)
    else:
        arg_parser = Tissue_expression_analyser(argparse_analysis())
    #This if statement is just to prevent any errors if no target tissues have been supplied
    if len(arg_parser.tissue) > 0:
        output_dir = arg_parser.outputdir
        #This function needs to run just for errorhandling, but the value is not needed and hence not stored
        arg_parser.check_meta_matrix_samples()
        #Initiating these object as None here allows use of if/else statements later on
        combined_specificity = None
        combined_enrichment = None
        #Looping over each target tissue allows individual specificity and enrichment analyses
        for tissue in arg_parser.tissue:
            specificity_result = arg_parser.run_specificity(tissue)
            enrichment_result = arg_parser.run_enrichment(tissue)
            #This if statement creates a file for each target tissue with boolean values of specificity and enrichment for each gene
            if specificity_result is not None or enrichment_result is not None:
                merged = pd.concat([specificity_result, enrichment_result], axis=1)
                merged.to_csv(f"{output_dir}/{tissue}_expression.tsv", sep='\t')

            #This if statement is used to store the combined boolean specificity values for all target tissues (essentially appending boolean values at each loop)
            if combined_specificity is None:
                combined_specificity = specificity_result
            else:
                #Using the bitwise operator "&" on boolean values returns True if both objects are True, but returns False if one or both of the objects are False
                #This way only genes that are specific to all target tissues remain True
                combined_specificity = combined_specificity & specificity_result
            
            #This if statement is used to store the combined boolean enrichment values for all target tissues (essentially appending boolean values at each loop)
            if combined_enrichment is None:
                combined_enrichment = enrichment_result
            else:
                #Using the bitwise operator "&" on boolean values returns True if both objects are True, but returns False if one or both of the objects are False
                #This way only genes that are enriched in all target tissues remain True
                combined_enrichment = combined_enrichment & enrichment_result

        #This if statement concatenates the pandas series of combined boolean specificity and enrichment into a pandas dataframe with a column for each and exporting this as a tsv file
        if specificity_result is not None or enrichment_result is not None:
            combined_merged = pd.concat([combined_specificity, combined_enrichment], axis=1)
            combined_merged.to_csv(f"{output_dir}/target_tissues_expression.tsv", sep='\t')
        
        #This if statement creates a tsv file with the genes found to be specific to all the target tissues
        if combined_specificity is not None:
            combined_specificity = combined_specificity[combined_specificity == True].index
            with open(f"{output_dir}/target_tissues_specificity.tsv", "w") as comb_spec_file:
                for gene in combined_specificity:
                    comb_spec_file.write(f"{gene}\n")
        
        #This if statement creates a tsv file with the genes found to be enriched in all the target tissues
        if combined_enrichment is not None:
            combined_enrichment = combined_enrichment[combined_enrichment == True].index
            with open(f"{output_dir}/target_tissues_enrichment.tsv", "w") as comb_enrich_file:
                for gene in combined_enrichment:
                    comb_enrich_file.write(f"{gene}\n")

    #The users input options are printet at the end of the info file for the user to refer back to if needed
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

    #The output directory is returned for when the script from another main script, for this subdirectory to be added to the main output directory
    return output_dir

#This bit makes it easier to run the analysis locally (useful for testing)
if __name__ == '__main__':
    execute_analysis()