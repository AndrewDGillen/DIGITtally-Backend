#Code written by Dr. Andrew Gillen (Dow-Davies lab, University of Glasgow)
#Originally written October 2022, last updated 26/06/2023

import os
import time

#Part of the pipeline for DIGITtally analysis - see www.digittally.org

#This program co-ordinates running ALL DIGITtally subordinate programs, given a "Context" dictionary, and the output folder name
#Please check subordinate programs for specifics/intricacies - Backbone only controls when/how these programs fire
#Descriptions of the ACTUAL PARAMETERS needed for each program are outlined in those programs 

#Parses score weighting options
#These are supplied in the format Metric:Weight, separated by ';'
def generic_score_weighting(score_string, source_context):

    for item in source_context:
        if "metricweight" in item:

            actual_metric = item.split('_')[1]
            metric_weight = source_context[item]

            score_string = score_string + f';{actual_metric}:{metric_weight}'
    
    return score_string

#Orthology weightstrings are constructed somewhat differently - each full metric has three subweights:
    #Enrichment (over whole insect)
    #Specificity (over non target tissues)
    #Any expression in target tissues
def build_orthology_weightstring(source_context):

    score_weightstring = f'Enrichment:{source_context["metricweight_Enrichment"]};Specificity:{source_context["metricweight_Specificity"]};Expression:{source_context["metricweight_Any"]}'
    
    return score_weightstring

#Formats data appropriately for FlyAtlas1 analysis using the FlyAtlas1_Analyser program
#n.b - when I refer to FlyAtlas1 I mean the original FlyAtlas (http://flyatlas.org/atlas.cgi)

#FlyAtlas1 Results will be stored in the FlyAtlas1_subfolder, demarkated with the supplied enrichment/abundance thresholds
def FA1(fa1_context, score_string, outputfolder):

    import FlyAtlas1_Analyser

    #Defines all variables from FlyAtlas1 based on user-defined context
    annotationfile = "AssociatedFiles/Affymetrix_Annotations.csv" 
    targettissues = fa1_context['tissues'] 

    try:
        permissive = fa1_context['permitted']
    except:
        permissive = []

    enrstringency = fa1_context['metricthreshold_Enrichment']
    abnstringency = fa1_context['metricthreshold_Specificity']
    mode = 'D'
    presentthresh = fa1_context['microarray_precenceThreshold']
    fa1file = 'AssociatedFiles/FlyAtlas1_Data.tsv'

    #The output folder for FlyAtlas1 output files is initialised
    outfolder = f'{outputfolder}/FlyAtlas1_Analysis_enr-' + str(enrstringency) + '_abn-' + str(abnstringency)
    os.mkdir(outfolder)

    #Runs the FlyAtlas1_Analyser. 
    #Enrichment/Abundance results are stored within the designated output folder. 
    FlyAtlas1_Analyser.analyse_FA1(annotationfile, targettissues, permissive, outfolder, enrstringency, abnstringency, mode, presentthresh, fa1file)

    fa1_location = outfolder

    #Next, we create a weight string determining how these results contribute to the FINAL DIGITtally
    #As with all data sources, this determines how INDIVIDUAL CHARACTERISTICS (eg Fly Type) contribute to final scoring
    #If a flytype is "OBLIGATORY"; ie the user has indicated a gene MUST pass thresholds in this type, that is marked imn the weightstring
    weightstring = 'MALE:1;FEMALE:1;'

    if fa1_context['type_use_Larval'] == True:
        larvalinuse = 1
    else:
        larvalinuse = 0

    weightstring = weightstring + f'LARVAL:{larvalinuse}'

    if fa1_context['type_obligatory_Larval'] == True:
        weightstring= weightstring + ':ob'
    
    if fa1_context['type_use_Adult'] == True:
        adultinuse = 1
    else:
        adultinuse = 0

    weightstring= weightstring + f';ADULTS:{adultinuse}'

    if fa1_context['type_obligatory_Adult'] == True:
        weightstring = weightstring +':ob'

    weightstring= weightstring + ';ALL:0'

    #If the user has indicated that they wish to use an extra threshold, that is enacted here
    if fa1_context["Enrichment_additionalThreshold"] == True or fa1_context["Specificity_additionalThreshold"] == True:

        extracheck = True
        print('processing additional thresholds')

        if fa1_context["Enrichment_additionalThreshold"] == True:
            enrstringency = fa1_context["Enrichment_additionalThresholdvalue"]
        
        if fa1_context["Specificity_additionalThreshold"] == True:
            abnstringency = fa1_context["Specificity_additionalThresholdvalue"]
        
        #Additional threshold data is held in its own output folder (Format same as previously)
        extrathresh_outfolder = f'{outputfolder}/FlyAtlas1_Analysis_extrathreshold_enr-' + str(enrstringency) + '_abn-' + str(abnstringency)

        os.mkdir(extrathresh_outfolder)

        FlyAtlas1_Analyser.analyse_FA1(annotationfile, targettissues, permissive, extrathresh_outfolder, enrstringency, abnstringency, mode, presentthresh, fa1file)

        fa1_extra_location = extrathresh_outfolder

    else:
        fa1_extra_location = ''
        extracheck = False
    
    #score weight dictionary is generated from the score_string
    score_weights = generic_score_weighting(score_string, fa1_context)

    #Output locations, weightstrings, score_weights are returned
    return fa1_location, fa1_extra_location, weightstring, score_weights, extracheck

#Formats data appropriately for FlyAtlas2 analysis using the FlyAtlas1_Analyser program
# FlyAtlas2 available at https://motif.mvls.gla.ac.uk/FlyAtlas2/

#FlyAtlas2 Results will be stored in the FlyAtlas2_subfolder, demarkated with the supplied enrichment/abundance thresholds
def FA2(fa2_context, sql_settings, score_string, outputfolder):

    import FlyAtlas2_Analyser

    #Defines all variables from FlyAtlas2 based on user-defined context
    database = 'DIGITtally_FA2_db'
    host = sql_settings['host']
    pw = sql_settings['pass']
    user = sql_settings['user']
    enrstringency = fa2_context['metricthreshold_Enrichment']
    abnstringency = fa2_context['metricthreshold_Specificity']
    weights = "FA2_Weighting"
    weightingmode = 1
    targettissues = fa2_context['tissues']

    try:
        permissive = fa2_context['permitted']
    except:
        permissive = []

    fpkmthresh = fa2_context['rnaseq_FPKMThreshold']

    #The output folder is initialised
    outfolder = f'{outputfolder}/FlyAtlas2_Analysis_enr-' + str(enrstringency) + '_abn-' + str(abnstringency)
    os.mkdir(outfolder)

    FlyAtlas2_Analyser.analyse_FA2(database, host, pw, user, enrstringency, abnstringency, weights, weightingmode, targettissues, permissive, fpkmthresh, outfolder)

    fa2_location = outfolder

    #Next, we create a weight string determining how these results contribute to the FINAL DIGITtally
    #As with all data sources, this determines how INDIVIDUAL CHARACTERISTICS (eg Fly Type) contribute to final scoring
    #If a flytype is "OBLIGATORY"; ie the user has indicated a gene MUST pass thresholds in this type, that is marked imn the weightstring
    weightstring = ''

    if fa2_context['type_use_Male'] == True:
        maleinuse = 1
    else:
        maleinuse = 0

    weightstring= weightstring + f'MALE:{maleinuse}'

    if fa2_context['type_obligatory_Male'] == True:
        weightstring = weightstring +':ob'
    
    if fa2_context['type_use_Female'] == True:
        femaleinuse = 1
    else:
        femaleinuse = 0

    weightstring = weightstring + f';FEMALE:{femaleinuse}'

    if fa2_context['type_obligatory_Female'] == True:
        weightstring = weightstring +':ob'

    if fa2_context['type_use_Larval'] == True:
        larvalinuse = 1
    else:
        larvalinuse = 0

    weightstring = weightstring + f';LARVAL:{larvalinuse}'

    if fa2_context['type_obligatory_Larval'] == True:
        weightstring = weightstring + ':ob'

    weightstring = weightstring + ';ADULTS:0;ALL:0'

    #If the user has indicated that they wish to use an extra threshold, that is enacted here
    if fa2_context["Enrichment_additionalThreshold"] == True or fa2_context["Specificity_additionalThreshold"] == True:
        print('processing additional thresholds')
        extracheck = True

        if fa2_context["Enrichment_additionalThreshold"] == True:
            enrstringency = fa2_context["Enrichment_additionalThresholdvalue"]
        
        if fa2_context["Specificity_additionalThreshold"] == True:
            abnstringency = fa2_context["Specificity_additionalThresholdvalue"]
        
        #Additional threshold data is held in its own output folder (Format same as previously)
        extrathresh_outfolder = f'{outputfolder}/FlyAtlas2_Analysis_extrathreshold_enr-' + str(enrstringency) + '_abn-' + str(abnstringency)

        os.mkdir(extrathresh_outfolder)

        FlyAtlas2_Analyser.analyse_FA2(database, host, pw, user, enrstringency, abnstringency, weights, weightingmode, targettissues, permissive, fpkmthresh, extrathresh_outfolder)

        fa2_extra_location = extrathresh_outfolder
    
    else:
        fa2_extra_location = ''
        extracheck = False
    
    #score weight dictionary is generated from the score_string
    score_weights = generic_score_weighting(score_string, fa2_context)

    #Output locations, weightstrings, score_weights are returned
    return fa2_location, fa2_extra_location, weightstring, score_weights, extracheck

#Uses FlyAtlas1 and FlyAtlas2 data to initialise a Tally output and gather a list of genes of interest
#This is ABSOLUTELY NECESSARY for handling FlyCellAtlas and FlyBase datasets, as these contain too much data to analyse for EVERY Drosophila gene
def StartTallying(context, fa1_location, fa1_extra_location, fa1weightstring, fa2_location, fa2_extra_location, fa2weightstring, outputfolder):

    import GOIAgglomerator

    #Defines all variables based on user-defined context

    #In particular, if the user has defined a list of genes they're interested in, we override automatic genelist ID with their list
    #As such, these genes will be scored across ALL metrics regardless of threshold passing
    if 'GeneList' in context:
        preexisting_list = context['GeneList']
    else:
        preexisting_list = []
    print(preexisting_list)

    fa1_folder = fa1_location
    fa2_folder = fa2_location
    fa1extra = fa1_extra_location
    fa2extra = fa2_extra_location

    #Defines which scoring metrics will be used for Genelist identification
    if 'Enrichment' in context['Scoring']['Metrics'] and 'Specificity' in context['Scoring']['Metrics']:
        mode = 'both'
    elif 'Enrichment' in context['Scoring']['Metrics']:
        mode = 'ENRICHED'
    elif 'Specificity' in context['Scoring']['Metrics']:
        mode = 'ABUNDANT'
    else:
        mode = "neither"


    #Initialises output folders for Genelists, their synonyms, and final Tally files
    os.mkdir(f"{outputfolder}/Genelists")
    os.mkdir(f"{outputfolder}/Genelists/Synonyms")
    os.mkdir(f'{outputfolder}/Tally_files')

    outputlists = f"{outputfolder}/Genelists"
    outputtally = f"{outputfolder}/Tally_files"

    infofile = "AssociatedFiles/FlyBase_GenesvTranscriptsvSymbols.tsv"
    secondaryannotations = "AssociatedFiles/FlyBase_SecondaryAnnotations.tsv"

    weightsfa1 = fa1weightstring
    weightsfa2 = fa2weightstring

    #Defines the maximum length of Gene of Interest lists - capped at 250 for feasibility in single cell analysis
    maxlen = 250
    
    #Initialises the tally file
    try:
        GOIAgglomerator.make_tally_basics(preexisting_list, fa1_folder, fa2_folder, fa1extra, fa2extra, mode, outputlists, outputtally, infofile, secondaryannotations, weightsfa1, weightsfa2, maxlen)
    except Exception as e:
        print(e)

#Checks whether the genes of interest, or their synonyms, can be identified in FlyCellAtlas data
#This is more efficient than determining this information at run-time during FlyCellAtlas analysis proper
def CheckSCOPE(reference_genes, outfolder, loomfile):
    
    import SCOPEscope

    #Defines all variables based on user-defined context
    genes_of_interest = f'{outfolder}/Genelists/DetectedGOIs_Symbol.txt'
    references = reference_genes
    output_folder = f'{outfolder}/Genelists'
    synonymfile = 'AssociatedFiles/FlyBase_Synonyms.tsv'

    SCOPEscope.CheckTheList(loomfile, genes_of_interest, references, output_folder, synonymfile)

#Determines gene of interest expression in each cell/tissue/annotation type within the FlyCellAtlas dataset
def FCA_Analysis_Proper(reference_gene_file, outfolder, loomfile):

    import FlyCellAtlas_Analyser

    #Defines all variables based on user-defined context
    goi_checked = f'{outfolder}/Genelists/DetectedGOIs_Symbol_checked.txt'
    reference_gene_file = reference_gene_file
    output_folder = f'{outfolder}/FlyCellAtlas_Analysis'

    #Initialises FCA output files
    os.mkdir(output_folder)
    os.mkdir(f'{output_folder}/OverlapPercents')
    os.mkdir(f'{output_folder}/OverlapCounts')

    FlyCellAtlas_Analyser.Analyse_FCA(loomfile, goi_checked, reference_gene_file, output_folder)
    
    return output_folder

#Uses FCA gene expression data to score each gene for Enrichment, Specificity, Ubiquity and Co-Expression (as selected by the user)
def FCA_TallyBuilder(celltypes, referencemode, outputfolder, fca_file_folder):

    import FlyCellAtlas_ToTally

    #Defines all variables based on user-defined context
    celltypes = celltypes
    synonyms = f'{outputfolder}/Genelists/DetectedGOIs_SynonymsForSCOPE.csv'
    outputatlas = f'{outputfolder}/FlyCellAtlas_Analysis'
    outputtally = f'{outputfolder}/Tally_files'
    referencemode = referencemode
    associatedfiles = 'AssociatedFiles'

    FlyCellAtlas_ToTally.build_the_tally(fca_file_folder, celltypes, synonyms, outputatlas, outputtally, referencemode, associatedfiles)

#High level function for co-ordinating FlyCellAtlas (FCA) analysis, the most complex analysis within DIGITtally
#The FlyCellAtlas project can be found at https://flycellatlas.org/
def FCA(fca_context, references, score_string, outputfolder,loomfile):

    #Checks whether genes of interest can be found within the FCA dataset
    CheckSCOPE(references, outputfolder, loomfile)
    
    print(references)
    
    #If the user has defined reference genes, we store those and turn on the reference mode
    if references != []:
        ref_file = f'{outputfolder}/Genelists/UserRefGenes_Symbol_checked.txt'
        run_ref = True
    else:
        ref_file = ''
        run_ref = False
    
    #On the DIGITtally site, FlyCellAtlas cell types are formatted for user readability. Here, we convert them back to machhine-readable format
    corrected_celltypes = []

    for celltype in fca_context['celltypes']:
        if '"' in celltype:
            celltype = celltype.replace('"', '')
        if "\\" in celltype:
            celltype = celltype.replace("\\", '')
        if "'" in celltype:
            celltype = celltype.replace("'", "")
        corrected_celltypes.append(celltype)
    
    #Gathers FCA expression data
    fca_folder_location = FCA_Analysis_Proper(ref_file, outputfolder, loomfile)
    
    #Scores genes based on FCA data
    FCA_TallyBuilder(corrected_celltypes, run_ref, outputfolder, fca_folder_location)
    
    #We return a weightstring informing the final analysis of how FCA data should be scored
    return generic_score_weighting(score_string, fca_context)

#Handles automated searching of the data repository FlyBase (www.flybase.org)
#FlyBase API is insufficient for this purpose (as it does not handle Controlled Vocabulary searching) so a Selenium Webdriver is utilised
#INTERNET ACCESS IS REQUIRED
def FB(fb_context, score_string, outputfolder):

    import FlyBase_Checker_neoCHROME

    #Defines all variables based on user-defined context
    input_fbids = f'{outputfolder}/Genelists/DetectedGOIs_FbID.txt'
    target_annotations = fb_context['vocab']
    outputfb = f'{outputfolder}/FlyBase_Analysis'
    outputtally = f'{outputfolder}/Tally_files'
    alleles_file = 'AssociatedFiles/FlyBase_Alleles.tsv'
    
    complete = 0
    errorcount = 0
    
    #We try up to 5 times to carry out a FlyBase search, to account for intermittent internet issues.
    while complete == 0 and errorcount < 5:
        try:
            FlyBase_Checker_neoCHROME.check_the_literature(input_fbids, target_annotations, outputfb, outputtally, alleles_file)
            complete = 1
        except Exception as e:
            time.sleep(60)
            print(e)
            errorcount += 1

    if errorcount >= 5:
        print('FLYBASE ERROR LIMIT REACHED')
        exit()
    return generic_score_weighting(score_string, fb_context)

#Searches DIOPT data for existance/disease association of Homo sapiens orthologs of Drosophila melanogaster genes of interest
#The DIOPT resource is available at https://www.flyrnai.org/cgi-bin/DRSC_orthologs.pl
def HSAP(hsap_context, outfolder):

    import Orthology_Hsap

    #Defines all variables based on user-defined context

    #In particular we detect which scoring metrics will be used
    input_fbids = f'{outfolder}/Genelists/DetectedGOIs_FbID.txt'

    if hsap_context['choice_UseOrthoDetection'] == True:
        detection_weight = hsap_context['metricweight_OrthoDetection']
    else:
        detection_weight = 0

    if hsap_context['choice_UseDiseaseAssociation'] == True:
        disease_weight = hsap_context['metricweight_DiseaseAssociation']
    else:
        disease_weight = 0

    stringency = hsap_context['DIOPT_orthoscoreThreshold']
    outputorthos = f'{outfolder}/Orthology'
    outputtally = f'{outfolder}/Tally_files'
    flybase_file = 'AssociatedFiles/FlyBase_HumanOrthology.tsv'

    #Carries out DIOPT analysis
    Orthology_Hsap.AnalyseDIOPTFindings(input_fbids, detection_weight, disease_weight, stringency, outputorthos, outputtally, flybase_file)

#Uses a flat file generated from the OrthoDB dataset to find genes of interest orthologs within other DIGITtally species
#OrthoDB is avaliable at https://www.orthodb.org/
def ORTHODB(species_list, outfolder):

    import OrthoDB_FromFile

    #Defines all variables based on user-defined context
    input_symbols = f'{outfolder}/Genelists/DetectedGOIs_FbID.txt'
    species = species_list
    output_folder = f'{outfolder}/Orthology'
    ortholog_file =  'AssociatedFiles/OrthoDB_Database_261022.csv'
    
    #Finds the OrthoDB orthologs
    OrthoDB_FromFile.gather_orthologs(input_symbols, ortholog_file, species, output_folder)

#Analyses MozAtlas (http://mozatlas.gen.cam.ac.uk/mozatlas/) and MozTubules (http://moztubules.org) Anopheles gambiae data 
#Uses OrthoDB-identified orthologs for the Gene of Interest list
def AGAM(mozatlas_context, moztubules_context, sql_settings,outputfolder):

    import Orthology_Agam

    #Defines all variables based on user-defined context
    database = 'DIGITtally_MozAtlas_db'
    host = sql_settings['host']
    pw = sql_settings['pass']
    user = sql_settings['user']

    orthologs = f'{outputfolder}/Orthology/Anopheles_gambiae_orthologs.csv'
    
    outputorthology = f'{outputfolder}/Orthology'
    outputtally = f'{outputfolder}/Tally_files'
    moztubules_probes = 'AssociatedFiles/MozTubules_Probe_IDs.csv'
    moztubules_expression = 'AssociatedFiles/MozTubules_ExpressionData.txt'
    
    modules = []

    #If the user wants to use MozAtlas, their settings are gathered and a score weight string is constructed
    if mozatlas_context != []:
        modules.append('MozAtlas')
        tissues = mozatlas_context['tissues']
        atlasthresh_enr =  mozatlas_context['metricthreshold_Enrichment']
        atlasthresh_spec = mozatlas_context['metricthreshold_Specificity']
        present_thresh = mozatlas_context['microarray_precenceThreshold']

        targets_weightstring = ''

        if mozatlas_context['type_use_Male'] == True:
            maleinuse = 1
        else:
            maleinuse = 0

        targets_weightstring = targets_weightstring + f'MALE:{maleinuse}'
        
        if mozatlas_context['type_use_Female'] == True:
            femaleinuse = 1
        else:
            femaleinuse = 0

        targets_weightstring = targets_weightstring + f';FEMALE:{femaleinuse};LARVAL:0;ADULT:0;ALL:0'

        ma_targetweights = targets_weightstring

        #Constructs a weight string which is used to adjust scores for the source based on user-selected importance of orthology measures
        ma_measureweights = build_orthology_weightstring(mozatlas_context)

    else:
        tissues = []
        atlasthresh_enr = 0
        atlasthresh_spec = 0
        present_thresh = 0

        ma_targetweights = 'MALE:0;FEMALE:0;LARVAL:0;ADULT:0;ALL:0'
        ma_measureweights = 'Enrichment:0;Abundance:0;Expression:0'

    #If the user wants to use MozTubules, their settings are gathered and a score weight string is constructed
    if moztubules_context != []:

        modules.append('MozTubules')
        tubulesthresh = moztubules_context['metricthreshold_Enrichment']
        tub_to_atlas = moztubules_context['moztubules_to_mozatlastissue']

        targets_weightstring = ''

        if moztubules_context['type_use_Male'] == True:
            maleinuse = 1
        else:
            maleinuse = 0

        targets_weightstring = targets_weightstring + f'MALE:{maleinuse}'
        
        if moztubules_context['type_use_Female'] == True:
            femaleinuse = 1
        else:
            femaleinuse = 0

        targets_weightstring = targets_weightstring + f';FEMALE:{femaleinuse}'

        if moztubules_context['type_use_Larval'] == True:
            larvalinuse = 1
        else:
            larvalinuse = 0

        targets_weightstring = targets_weightstring + f';LARVAL:{larvalinuse}'
        
        if moztubules_context['type_use_Adult'] == True:
            adultinuse = 1
        else:
            adultinuse = 0

        targets_weightstring = targets_weightstring + f';ADULT:{adultinuse};ALL:0'

        mt_targetweights = targets_weightstring
        
    else:
        tubulesthresh = 0
        tub_to_atlas = 0
        
        mt_targetweights = 'MALE:0;FEMALE:0;LARVAL:0;ADULT:0;ALL:0'
    
    #The Orthology_Agam module is called which analyses ALL Anopheles gambiae data (Ie both MozTubules and MozAtlas) 
    Orthology_Agam.Score_Agam_Orthology(
        database, host, pw, user, orthologs, outputorthology, outputtally, moztubules_probes, moztubules_expression, tissues, atlasthresh_enr, 
        atlasthresh_spec, present_thresh, ma_targetweights, ma_measureweights, tubulesthresh, tub_to_atlas, mt_targetweights, modules
        )

#Analyses Aegypti-Atlas (http://aegyptiatlas.buchonlab.com/) Aedes aegypti data 
#Uses OrthoDB-identified orthologs for the Gene of Interest list
def AAEG(aegypti_context, outputfolder):

    import Orthology_Aaeg

    #Defines all variables based on user-defined context
    orthologs = f'{outputfolder}/Orthology/Aedes_aegypti_orthologs.csv'
    
    enrichthresh = aegypti_context['metricthreshold_Enrichment']
    specthresh = aegypti_context['metricthreshold_Specificity']
    tissues = aegypti_context['tissues']
    aegyptiatlas_file = 'AssociatedFiles/AegyptiAtlas_FPKMs.csv'
    outputorthology = f'{outputfolder}/Orthology'
    outputtally = f'{outputfolder}/Tally_files'
    fpkm_background = aegypti_context['rnaseq_FPKMThreshold']

    #Constructs a weight string which is used to adjust scores for the source based on user-selected importance of orthology measures
    measureweights = build_orthology_weightstring(aegypti_context)

    #Analyses A.aeg expression data
    Orthology_Aaeg.Score_Aaeg_Orthology(orthologs, enrichthresh, specthresh, tissues, aegyptiatlas_file, outputorthology, outputtally, fpkm_background, measureweights)

#SilkDB contains many Silkworm life stages.
#For digittally.org, these life stages need to be presented in a user-readable format
#Here, they are converted back into a machine-readable format
def handle_bmor_life_stages(context):

    stage_style_conversion = {
        'Larvae_4thInstar' : '4th-instar-day-3',
        'Larvae_4thMolt' : 'fourth-larval-molting',
        'Larvae_5thInstar':['5th-instar-day-0:1', '5th-instar-day-3'],
        'Larvae_Wandering':'wandering',
        'Larvae_PrePupa' : 'pre-pupa',
        'Pupa': ['pupa-day-1', 'pupa-day-4', 'pupa-day-7-8'],
        'Moth':'moth-day-1'
    }

    bmor_string = ''

    for subcontext in context:
        if'type_use_' in subcontext:
            silkwormtype = subcontext[9:]
            
            if context[subcontext] == True:
                use = 1
            else:
                use = 0
            
            if type(stage_style_conversion[silkwormtype]) == list:
                for correspondingtype in stage_style_conversion[silkwormtype]:
                    bmor_string = bmor_string + f'{correspondingtype}:{use};'
            else:
                bmor_string = bmor_string + f'{stage_style_conversion[silkwormtype]}:{use};'
    
    bmor_string = bmor_string[:-1]

    return bmor_string

#Analyses SilkDB (https://silkdb.bioinfotoolkits.net/) Bombyx mori data 
#Uses OrthoDB-identified orthologs for the Gene of Interest list
def BMOR(bmor_context, outputfolder):

    import Orthology_Bmor

    #Defines all variables based on user-defined context
    ortholog_file = f'{outputfolder}/Orthology/Bombyx_mori_orthologs.csv'
    
    enrichthresh = bmor_context['metricthreshold_Enrichment']
    specthresh = bmor_context['metricthreshold_Specificity']
    tissues = bmor_context['tissues']
    silkdb_file = 'AssociatedFiles/Silkdb_Expression.txt'
    outputorthology = f'{outputfolder}/Orthology'
    outputtally = f'{outputfolder}/Tally_files'
    fpkm_background = bmor_context['rnaseq_FPKMThreshold']

    #Gathers user-defined weighting for individual silkworm lifestages
    typeweights = handle_bmor_life_stages(bmor_context)

    #Sets a switch for the "Any Stage" approach
    #Ie if a gene passes Enrichmet/Abundance/Expression thresholds at ANY life stage, it earns a point
    if bmor_context['choice_UseAnyStage'] == True:
        stageweighting = 'Any_type:1;'
    else:
        stageweighting = 'Any_type:0;'
    
    #Sets a switch for the "ALL Stage" approach
    #Ie if a gene passes Enrichmet/Abundance/Expression thresholds at ALL life stages, it earns a point
    #NOTE: This can be used in conjunction with the previous setting
    if bmor_context['choice_UseAllStage'] == True:
        stageweighting = stageweighting + 'All_type:1;Individual_stage:0'
    else:
        stageweighting = stageweighting + 'All_type:0;Individual_stage:0'

    lifestageweight = stageweighting

    #Constructs a weight string which is used to adjust scores for the source based on user-selected importance of orthology measures
    measureweights = build_orthology_weightstring(bmor_context)

    #Analyses B.mor expression data
    Orthology_Bmor.Score_Bmor_orthology(ortholog_file, enrichthresh, specthresh, tissues, silkdb_file, outputorthology, outputtally, fpkm_background, typeweights, lifestageweight, measureweights)

#A simple function to return weights for a specific datasource
def return_global_weight(context, datasource):
    return context['Sources_global']['Weights'][f'{datasource}']

#Builds a "GLOBAL SCORE" weightstring, containing weightings (defined by user) for each specific score type
def gather_global_weights(context):
    broad_score_string = ''
    allmetrics = ['Enrichment','Specificity','Ubiquity','Coexpression','Known Activity','Orthology']

    for scoremetric in allmetrics:
        try:
            broad_score_string = broad_score_string + f'{scoremetric}:{context["Scoring"]["Weights"][scoremetric]};'
        except:
            broad_score_string = broad_score_string + f'{scoremetric}:0;'
    
    broad_score_string = broad_score_string[:-1]

    return broad_score_string

#Gathers:
#  all data sources in use, 
# their relative importance (based on user input), 
# The relative importance of each scoring metric GLOBALLY,
# The relative importance of each scoring metric within each data source

#And uses analysed data to construct the final "DIGITtally"
def CompleteAnalysis(global_score_weights, fa1weights, fa2weights, fcaweights, fbweights, orthoweights, orthospecies,outputfolder):

    import Finalise_DIGITtally

    #Gathers user defined weights
    tally_file_folder = f'{outputfolder}/Tally_files'
    synonym_file = f'{outputfolder}/Genelists/Synonyms/DetectedGOIs_Synonyms.csv'
    broad_categories = global_score_weights
    flyatlas1 = fa1weights
    flyatlas2 = fa2weights
    flycellatlas = fcaweights
    flybase = fbweights
    orthology = orthoweights
    species = orthospecies

    #Creates a DIGITtally for the user-defined settings
    Finalise_DIGITtally.WrapUp(tally_file_folder, synonym_file, broad_categories, flyatlas1, flyatlas2, flycellatlas, flybase, orthology, species)

#This function co-ordinates the parts of DIGITtally, running components as requested and collecting variables for the final tally calculations
#Note this program prints to terminal - this is useful for error logging
def run_from_backbone(context, folder_addendum):

    from selenium.common.exceptions import NoSuchElementException, TimeoutException
    import json
    import sys


    #Gathers information from the CONFIG file
    #FIle location is not shown here for data security reasons
    with open('CONFIG FILE') as config_file:
        config = json.load(config_file)

    #Run information (time taken, errors etc) is logged to a seperate LOGGING folder
    #This helps with accessing user information in the case of any issues
    with open(f'{config["LogLocation"]}/DT_Log_{folder_addendum}.txt', 'w+') as sys.stdout:
        print('RUNNING!')
        try:
            start_time = time.time()
            last_time = start_time
        except Exception as e:
            print(e)

        print(f'Beginning handling request {folder_addendum}')
        outputfolder = f"DIGITtally_results_{folder_addendum}"
        os.mkdir(outputfolder)
        
        #This is important to define here to prevent loompy from locking .loom files, ensuring parallel processes can be run
        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
        
        #The user's chosen sources are gathered into a list
        sources = context["Sources_global"]["Sources"]

        #Orthology settings are defined
        orthology_sources = ["MozAtlas", "MozTubules", "AegyptiAtlas", "SilkDB"]
        orthology_species = {
            'Hsap': 'Homo sapiens',
            'Agam':'Anopheles gambiae',
            'Aaeg':'Aedes aegypti', 
            'Bmor':'Bombyx mori'
        }

        #Generic settings for connecting to MySQL
        sql_settings = {
            'host':config['SQL_HOST'],
            'user':config['SQL_USER'],
            'pass':config['SQL_PASS']
        }

        extra_threshold_sources = []

        #Handles FlyAtlas1, passes needed variables to the Agglomerator
        if "FlyAtlas1" in sources:
            
            print('Processing FlyAtlas1!')
            fa1_metricweights = f'GLOBAL:{return_global_weight(context, "FlyAtlas1")}'
            fa1_location, fa1_extra_location, fa1weightstring, fa1_metricweights, extracheck = FA1(context['FlyAtlas1'], fa1_metricweights, outputfolder)

            if extracheck == True:
                extra_threshold_sources.append('FlyAtlas1')
            
            print('\n FA1 DONE - %s seconds' % (time.time() - last_time))
            last_time = time.time()

        else:
            fa1_location = ''
            fa1_extra_location = ''
            fa1weightstring = "MALE:0;FEMALE:0;lARVAL:0;ADULTS:0;ALL:0"
            fa1_metricweights = 'GLOBAL:0;Enrichment:0;Specificity:0'
        
        #Handles FlyAtlas2, passes needed variables to the Agglomerator
        if "FlyAtlas2" in sources:
            print('Processing FlyAtlas2!')
            fa2_metricweights = f'GLOBAL:{return_global_weight(context, "FlyAtlas2")}'
            try:
                fa2_location, fa2_extra_location, fa2weightstring, fa2_metricweights, extracheck = FA2(context['FlyAtlas2'], sql_settings, fa2_metricweights, outputfolder)
            except Exception as e:
                print(e)
                exit()
            if extracheck == True:
                extra_threshold_sources.append('FlyAtlas2')
            
            print('\n FA2 DONE - %s seconds' % (time.time() - last_time))
            last_time = time.time()

        else:
            fa2_location = ''
            fa2_extra_location = ''
            fa2weightstring = "MALE:0;FEMALE:0;lARVAL:0;ADULTS:0;ALL:0"
            fa2_metricweights = 'GLOBAL:0;Enrichment:0;Specificity:0'

        #Initialises the DIGITtally "TALLY" file
        #This is the basis for the final output, and is used to define genes of interest based on FlyAtlas1 and FlyAtlas2 data
        try:
            print('Initialising the Tally')
            StartTallying(context, fa1_location, fa1_extra_location, fa1weightstring, fa2_location, fa2_extra_location, fa2weightstring, outputfolder)
        except Exception as e:
            print(e)

        #Handles FlyCellAtlas, passes needed variables to the Agglomerator
        if "FlyCellAtlas" in sources:

            print('Processing FlyCellAtlas!')

            #Gathers user-defined Reference genes 
            if 'References' in context:
                ref_genes = context['References']
            else:
                ref_genes = ''

            #
            fca_metricweights = f'GLOBAL:{return_global_weight(context, "FlyCellAtlas")}'
            
            loom_inc = 0
            finished = 0

            import time

            #This try/except block iterates through .loom files (I keep ten copies within Associated Files)
            #If a given .loom can't be accessed, we wait a minute and try the next one
            while loom_inc < 10 and finished == 0:
                loomfile = f'AssociatedFiles/FCA_WholeFly_Loom_{loom_inc}.loom'

                try:
                    fca_metricweights = FCA(context['FlyCellAtlas'], ref_genes, fca_metricweights, outputfolder,loomfile)
                    finished = 1

                #This error can occur if a given LOOM file is locked.
                except BlockingIOError:
                    loom_inc += 1
                    print('ERROR HIT')

                    if loom_inc == 10:
                        time.sleep(60)
                        loom_inc = 0

            print('\n FCA DONE - %s seconds' % (time.time() - last_time))
            last_time = time.time()
        else:
            fca_metricweights = 'GLOBAL:0;Enrichment:0;Specificity:0;Ubiquity:0;Coexpression:0'

        #Handles FlyBase, passes needed variables to the Agglomerator
        if "FlyBase" in sources:
            print('Processing FlyBase!')
            
            #Initialises folders for processing FlyBase data
            os.mkdir(f'{outputfolder}/FlyBase_Analysis')
            os.mkdir(f'{outputfolder}/FlyBase_Analysis/Downloads')
            os.mkdir(f'{outputfolder}/FlyBase_Analysis/Downloaded_lists')
            
            done = 0
            
            #This can error out if an internet issue exists server-side. These errors are tracked.
            #Strictly, the try/except is not currently necessary - it used to be used to retry under certain circumstances
            while done == 0:
                try:
                    fb_metricweights = f'GLOBAL:{return_global_weight(context, "FlyBase")}'
            
                    fb_metricweights = FB(context['FlyBase'], fb_metricweights, outputfolder)
                    done = 1
                    
                except Exception as e:
                    print(e)
                    exit()   

            print('\n FLYBASE DONE - %s seconds' % (time.time() - last_time))
            last_time = time.time()

        else:
            fb_metricweights = 'GLOBAL:0;Association:0;Phenotype:0'

        species = []

        #Orthology is handled in its own block, as users may choose to exclude ALL orthology methods broadly
        if 'Orthology' in context['Scoring']['Metrics']:
            
            speciesweightstring = f'GLOBAL:{context["Scoring"]["Weights"]["Orthology"]}'

            #An Orthology folder is initialised
            os.mkdir(f'{outputfolder}/Orthology')
            
            tmpspecies = []
            

            for ortho_species in orthology_species:

                species_weight = context["Scoring"]["Weights"][ortho_species]
                speciesweightstring = speciesweightstring + f';{ortho_species}:{species_weight}'
                
                if context["Scoring"]["Weights"][ortho_species] != 0:
                    tmpspecies.append(orthology_species[ortho_species])
            
            print(tmpspecies)
            
            #Handles DIOPT, passes needed variables to the Agglomerator
            if "DIOPT" in sources:

                print('Processing DIOPT!')
                species.append('Homo sapiens')

                HSAP(context['DIOPT'], outputfolder)

                print('\n DIOPT DONE - %s seconds' % (time.time() - last_time))
                last_time = time.time()

            #Gathers Orthologs for each gene of interest from all species the user wishes to use
            if any(item in sources for item in orthology_sources):

                if 'Homo sapiens' in tmpspecies:
                    tmpspecies.remove('Homo sapiens')

                print('Finding Orthologs!')
                
                ORTHODB(tmpspecies, outputfolder)

                print('\n ORTHOLOG ID DONE - %s seconds' % (time.time() - last_time))
                last_time = time.time()

            #Handles MozAtlas, passes needed variables to the Agglomerator
            if "MozAtlas" in sources or "MozTubules" in sources:

                print('Processing Anopheles gambiae!')
                species.append('Anopheles gambiae')

                if "MozAtlas" in sources:
                    ma_context = context['MozAtlas']
                else:
                    ma_context = []

                if "MozTubules" in sources:
                    mt_context = context['MozTubules']
                else:
                    mt_context = []
                
                AGAM(ma_context, mt_context, sql_settings, outputfolder)

                print('\n AGAM DONE - %s seconds' % (time.time() - last_time))
                last_time = time.time()

            #Handles Aegypti-Atlas, passes needed variables to the Agglomerator
            if "AegyptiAtlas" in sources:

                print('Processing AegyptiAtlas!')
                species.append('Aedes aegypti')

                AAEG(context['AegyptiAtlas'], outputfolder)

                print('\n AAEG DONE - %s seconds' % (time.time() - last_time))
                last_time = time.time()

            #Handles SilkDB, passes needed variables to the Agglomerator
            if "SilkDB" in sources:

                print('Processing SilkDB!')
                species.append('Bombyx mori')

                BMOR(context['SilkDB'], outputfolder)

                print('\n BMORI DONE - %s seconds' % (time.time() - last_time))
                last_time = time.time()

        else:
            speciesweightstring = 'GLOBAL:0;Hsap:0;Agam:0;Aaeg:0;Bmor:0'

        #Weights for each score type are gathered
        global_score_weights = gather_global_weights(context)
        
        print('Finishing the Analysis!')
        print(fb_metricweights)

        #Summarises all data into a final score - the DIGITtally
        try:
            CompleteAnalysis(global_score_weights, fa1_metricweights, fa2_metricweights, fca_metricweights, fb_metricweights, speciesweightstring, species, outputfolder)
        except Exception as e:
            print(e)

        print('\n FINAL TALLYING DONE - %s seconds' % (time.time() - last_time))
        print('\n\n DIGITtally DONE - %s seconds TOTAL' % (time.time() - start_time))