import os


def generic_score_weighting(score_string, source_context):

    for item in source_context:
        if "metricweight" in item:

            actual_metric = item.split('_')[1].title()
            metric_weight = source_context[item]

            if metric_weight == None:
                metric_weight = 0
            score_string = score_string + f';{actual_metric}:{metric_weight}'
    
    return score_string

def build_orthology_weightstring(source_context):

    score_weightstring = f'Enrichment:{source_context["metricweight_Enrichment"]};Specificity:{source_context["metricweight_Specificity"]};Expression:{source_context["metricweight_Any"]}'
    
    return score_weightstring

def UserUpload(usr_context, usr_file_locations, score_string, outputfolder, config):

    try:
        import MS_UserSample_Analyse
    except Exception as e:
        print(e)
        exit(1)

    from argparse import Namespace

    if usr_context['flytypes_selected'] == ["NA"]:
        sex_submission = []
    else:
        sex_submission=usr_context['flytypes_selected']

    if usr_context['ages_selected'] == ["NA"]:
        age_submission = []
    else:
        age_submission=usr_context['ages_selected']

    print("Imports OK")
    #We manually define arguments for the UserSample analysis program
    arg_object = Namespace(
        metadata=usr_file_locations['meta_matrix'],
        matrix=usr_file_locations['expression_matrix'],
        s=usr_context['use_usr_specificity'],
        sthres=usr_context['metricthreshold_specificity'],
        e=usr_context['use_usr_enrichment'],
        ethres=usr_context['metricthreshold_enrichment'],
        bthres=usr_context['metricthreshold_background'],
        tissue=usr_context['tissues_selected'],
        permissive=usr_context['permitted'],
        directory=outputfolder,
        age=age_submission,
        sex=sex_submission,
        logdir=config["LogLocation"],
        decimal=usr_context['usr_decimalpoint']
        )
    
    print("starting analysis")

    try:
        usr_outloc = MS_UserSample_Analyse.execute_analysis(artificial_args=arg_object) 
        usr_score_weights = generic_score_weighting(score_string, usr_context)
    except Exception as e:
        print(e)
        exit(1)
    
    print(usr_score_weights)

    return usr_score_weights, usr_outloc

def FA1(fa1_context, score_string, outputfolder):

    import FlyAtlas1_Analyser

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

    outfolder = f'{outputfolder}/FlyAtlas1_Analysis_enr-' + str(enrstringency) + '_abn-' + str(abnstringency)
    os.mkdir(outfolder)

    FlyAtlas1_Analyser.analyse_FA1(annotationfile, targettissues, permissive, outfolder, enrstringency, abnstringency, mode, presentthresh, fa1file)

    fa1_location = outfolder

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

    if fa1_context["Enrichment_additionalThreshold"] == True or fa1_context["Specificity_additionalThreshold"] == True:

        extracheck = True
        print('processing additional thresholds')

        if fa1_context["Enrichment_additionalThreshold"] == True:
            enrstringency = fa1_context["Enrichment_additionalThresholdvalue"]
        
        if fa1_context["Specificity_additionalThreshold"] == True:
            abnstringency = fa1_context["Specificity_additionalThresholdvalue"]
        
        extrathresh_outfolder = f'{outputfolder}/FlyAtlas1_Analysis_extrathreshold_enr-' + str(enrstringency) + '_abn-' + str(abnstringency)

        os.mkdir(extrathresh_outfolder)

        FlyAtlas1_Analyser.analyse_FA1(annotationfile, targettissues, permissive, extrathresh_outfolder, enrstringency, abnstringency, mode, presentthresh, fa1file)

        fa1_extra_location = extrathresh_outfolder
    
    else:
        fa1_extra_location = ''
        extracheck = False
    
    score_weights = generic_score_weighting(score_string, fa1_context)

    return fa1_location, fa1_extra_location, weightstring, score_weights, extracheck

def FA2(fa2_context, sql_settings, score_string, outputfolder):

    import FlyAtlas2_Analyser

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

    outfolder = f'{outputfolder}/FlyAtlas2_Analysis_enr-' + str(enrstringency) + '_abn-' + str(abnstringency)

    os.mkdir(outfolder)

    FlyAtlas2_Analyser.analyse_FA2(database, host, pw, user, enrstringency, abnstringency, weights, weightingmode, targettissues, permissive, fpkmthresh, outfolder)

    fa2_location = outfolder

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

    if fa2_context["Enrichment_additionalThreshold"] == True or fa2_context["Specificity_additionalThreshold"] == True:
        print('processing additional thresholds')
        extracheck = True

        if fa2_context["Enrichment_additionalThreshold"] == True:
            enrstringency = fa2_context["Enrichment_additionalThresholdvalue"]
        
        if fa2_context["Specificity_additionalThreshold"] == True:
            abnstringency = fa2_context["Specificity_additionalThresholdvalue"]
        
        extrathresh_outfolder = f'{outputfolder}/FlyAtlas2_Analysis_extrathreshold_enr-' + str(enrstringency) + '_abn-' + str(abnstringency)

        os.mkdir(extrathresh_outfolder)

        FlyAtlas2_Analyser.analyse_FA2(database, host, pw, user, enrstringency, abnstringency, weights, weightingmode, targettissues, permissive, fpkmthresh, extrathresh_outfolder)

        fa2_extra_location = extrathresh_outfolder
    
    else:
        fa2_extra_location = ''
        extracheck = False
    
    score_weights = generic_score_weighting(score_string, fa2_context)

    return fa2_location, fa2_extra_location, weightstring, score_weights, extracheck

def StartTallying(context, fa1_location, fa1_extra_location, fa1weightstring, fa2_location, fa2_extra_location, fa2weightstring, outputfolder, usr_folder, usr_weights):
    import GOIAgglomerator

    if 'GeneList' in context:
        preexisting_list = context['GeneList']
    else:
        preexisting_list = []
    print(preexisting_list)
    fa1_folder = fa1_location
    fa2_folder = fa2_location
    fa1extra = fa1_extra_location
    fa2extra = fa2_extra_location

    if 'Enrichment' in context['Scoring']['Metrics'] and 'Specificity' in context['Scoring']['Metrics']:
        mode = 'both'
    elif 'Enrichment' in context['Scoring']['Metrics']:
        mode = 'ENRICHED'
    elif 'Specificity' in context['Scoring']['Metrics']:
        mode = 'ABUNDANT'
    else:
        mode = "neither"
    
    os.mkdir(f"{outputfolder}/Genelists")
    os.mkdir(f"{outputfolder}/Genelists/Synonyms")
    os.mkdir(f'{outputfolder}/Tally_files')

    outputlists = f"{outputfolder}/Genelists"
    outputtally = f"{outputfolder}/Tally_files"

    infofile = "AssociatedFiles/FlyBase_GenesvTranscriptsvSymbols.tsv"
    secondaryannotations = "AssociatedFiles/FlyBase_SecondaryAnnotations.tsv"
    synonymfile = 'AssociatedFiles/FlyBase_Synonyms.tsv'

    weightsfa1 = fa1weightstring
    weightsfa2 = fa2weightstring
    maxlen = 250
    
    try:
        GOIAgglomerator.make_tally_basics(
            preexisting_list, 
            fa1_folder, 
            fa2_folder, 
            usr_folder, 
            fa1extra, 
            fa2extra, 
            mode, 
            outputlists, 
            outputtally, 
            infofile, 
            secondaryannotations, 
            weightsfa1, 
            weightsfa2, 
            usr_weights, 
            maxlen,
            synonymfile,
            )
    except Exception as e:
        print(e)

def CheckSCOPE(reference_genes, outfolder, loomfile):
    
    import SCOPEscope

    genes_of_interest = f'{outfolder}/Genelists/DetectedGOIs_Symbol.txt'
    references = reference_genes
    output_folder = f'{outfolder}/Genelists'
    synonymfile = 'AssociatedFiles/FlyBase_Synonyms.tsv'

    SCOPEscope.CheckTheList(loomfile, genes_of_interest, references, output_folder, synonymfile)

def FCA_Analysis_Proper(reference_gene_file, outfolder, loomfile):
    import FlyCellAtlas_Analyser

    goi_checked = f'{outfolder}/Genelists/DetectedGOIs_Symbol_checked.txt'
    reference_gene_file = reference_gene_file
    output_folder = f'{outfolder}/FlyCellAtlas_Analysis'

    os.mkdir(output_folder)
    os.mkdir(f'{output_folder}/OverlapPercents')
    os.mkdir(f'{output_folder}/OverlapCounts')

    FlyCellAtlas_Analyser.Analyse_FCA(loomfile, goi_checked, reference_gene_file, output_folder)
    
    return output_folder

def FCA_TallyBuilder(celltypes, referencemode, outputfolder, fca_file_folder):
    import FlyCellAtlas_ToTally

    celltypes = celltypes
    synonyms = f'{outputfolder}/Genelists/DetectedGOIs_SynonymsForSCOPE.csv'
    outputatlas = f'{outputfolder}/FlyCellAtlas_Analysis'
    outputtally = f'{outputfolder}/Tally_files'
    referencemode = referencemode
    associatedfiles = 'AssociatedFiles'

    FlyCellAtlas_ToTally.build_the_tally(fca_file_folder, celltypes, synonyms, outputatlas, outputtally, referencemode, associatedfiles)

def FCA(fca_context, references, score_string, outputfolder,loomfile):
    CheckSCOPE(references, outputfolder, loomfile)
    
    print(references)
    
    if references != []:
        ref_file = f'{outputfolder}/Genelists/UserRefGenes_Symbol_checked.txt'
        run_ref = True
    else:
        ref_file = ''
        run_ref = False
    
    corrected_celltypes = []

    for celltype in fca_context['celltypes']:
        if '"' in celltype:
            celltype = celltype.replace('"', '')
        if "\\" in celltype:
            celltype = celltype.replace("\\", '')
        if "'" in celltype:
            celltype = celltype.replace("'", "")
        corrected_celltypes.append(celltype)
    
    fca_folder_location = FCA_Analysis_Proper(ref_file, outputfolder, loomfile)
    
    FCA_TallyBuilder(corrected_celltypes, run_ref, outputfolder, fca_folder_location)
    
    return generic_score_weighting(score_string, fca_context)

def FB(fb_context, score_string, outputfolder):

    import FlyBase_FireFoxUpdate

    input_fbids = f'{outputfolder}/Genelists/DetectedGOIs_FbID.txt'
    target_annotations = fb_context['vocab']
    outputfb = f'{outputfolder}/FlyBase_Analysis'
    outputtally = f'{outputfolder}/Tally_files'
    alleles_file = 'AssociatedFiles/FlyBase_Alleles.tsv'
    
    complete = 0
    errorcount = 0
    
    while complete == 0 and errorcount < 5:
        try:
            FlyBase_FireFoxUpdate.check_the_literature(input_fbids, target_annotations, outputfb, outputtally, alleles_file)
            complete = 1
        except Exception as e:
            print(e)
            errorcount += 1
    if errorcount >= 5:
        raise ValueError

    return generic_score_weighting(score_string, fb_context)

def HSAP(hsap_context, outfolder):

    import Orthology_Hsap

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

    Orthology_Hsap.AnalyseDIOPTFindings(input_fbids, detection_weight, disease_weight, stringency, outputorthos, outputtally, flybase_file)

def ORTHODB(species_list, outfolder):

    import OrthoDB_FromFile

    input_symbols = f'{outfolder}/Genelists/DetectedGOIs_FbID.txt'
    species = species_list
    output_folder = f'{outfolder}/Orthology'
    ortholog_file =  'AssociatedFiles/OrthoDB_Database_261022.csv'
    
    OrthoDB_FromFile.gather_orthologs(input_symbols, ortholog_file, species, output_folder)
    
def AGAM(mozatlas_context, moztubules_context, sql_settings,outputfolder):

    import Orthology_Agam

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

        ma_measureweights = build_orthology_weightstring(mozatlas_context)

    else:
        tissues = []
        atlasthresh_enr = 0
        atlasthresh_spec = 0
        present_thresh = 0

        ma_targetweights = 'MALE:0;FEMALE:0;LARVAL:0;ADULT:0;ALL:0'
        ma_measureweights = 'Enrichment:0;Abundance:0;Expression:0'

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
        
    Orthology_Agam.Score_Agam_Orthology(
        database, host, pw, user, orthologs, outputorthology, outputtally, moztubules_probes, moztubules_expression, tissues, atlasthresh_enr, 
        atlasthresh_spec, present_thresh, ma_targetweights, ma_measureweights, tubulesthresh, tub_to_atlas, mt_targetweights, modules
        )

def AAEG(aegypti_context, outputfolder):

    import Orthology_Aaeg

    orthologs = f'{outputfolder}/Orthology/Aedes_aegypti_orthologs.csv'
    
    enrichthresh = aegypti_context['metricthreshold_Enrichment']
    specthresh = aegypti_context['metricthreshold_Specificity']
    tissues = aegypti_context['tissues']
    aegyptiatlas_file = 'AssociatedFiles/AegyptiAtlas_FPKMs.csv'
    outputorthology = f'{outputfolder}/Orthology'
    outputtally = f'{outputfolder}/Tally_files'
    fpkm_background = aegypti_context['rnaseq_FPKMThreshold']

    measureweights = build_orthology_weightstring(aegypti_context)

    Orthology_Aaeg.Score_Aaeg_Orthology(orthologs, enrichthresh, specthresh, tissues, aegyptiatlas_file, outputorthology, outputtally, fpkm_background, measureweights)

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

    string = ''

    for subcontext in context:
        if'type_use_' in subcontext:
            silkwormtype = subcontext[9:]
            
            if context[subcontext] == True:
                use = 1
            else:
                use = 0
            
            if type(stage_style_conversion[silkwormtype]) == list:
                for correspondingtype in stage_style_conversion[silkwormtype]:
                    string = string + f'{correspondingtype}:{use};'
            else:
                string = string + f'{stage_style_conversion[silkwormtype]}:{use};'
    
    string = string[:-1]

    return string

def BMOR(bmor_context, outputfolder):

    import Orthology_Bmor

    ortholog_file = f'{outputfolder}/Orthology/Bombyx_mori_orthologs.csv'
    
    enrichthresh = bmor_context['metricthreshold_Enrichment']
    specthresh = bmor_context['metricthreshold_Specificity']
    tissues = bmor_context['tissues']
    silkdb_file = 'AssociatedFiles/Silkdb_Expression.txt'
    outputorthology = f'{outputfolder}/Orthology'
    outputtally = f'{outputfolder}/Tally_files'
    fpkm_background = bmor_context['rnaseq_FPKMThreshold']

    typeweights = handle_bmor_life_stages(bmor_context)

    if bmor_context['choice_UseAnyStage'] == True:
        stageweighting = 'Any_type:1;'
    else:
        stageweighting = 'Any_type:0;'
    
    if bmor_context['choice_UseAllStage'] == True:
        stageweighting = stageweighting + 'All_type:1;Individual_stage:0'
    else:
        stageweighting = stageweighting + 'All_type:0;Individual_stage:0'

    lifestageweight = stageweighting

    measureweights = build_orthology_weightstring(bmor_context)

    Orthology_Bmor.Score_Bmor_orthology(ortholog_file, enrichthresh, specthresh, tissues, silkdb_file, outputorthology, outputtally, fpkm_background, typeweights, lifestageweight, measureweights)

def return_global_weight(context, datasource):
    return context['Sources_global']['Weights'][f'{datasource}']

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

def CompleteAnalysis(global_score_weights, usr_weights, fa1weights, fa2weights, fcaweights, fbweights, orthoweights, orthospecies,outputfolder):
    import Finalise_DIGITtally

    tally_file_folder = f'{outputfolder}/Tally_files'
    synonym_file = f'{outputfolder}/Genelists/Synonyms/DetectedGOIs_Synonyms.csv'
    broad_categories = global_score_weights
    flyatlas1 = fa1weights
    flyatlas2 = fa2weights
    flycellatlas = fcaweights
    flybase = fbweights
    orthology = orthoweights
    species = orthospecies

    Finalise_DIGITtally.WrapUp(tally_file_folder, synonym_file, broad_categories, usr_weights, flyatlas1, flyatlas2, flycellatlas, flybase, orthology, species)

#This function co-ordinates the parts of DIGITtally, running components as requested and collecting variables for the final tally calculations
def run_from_backbone(context, folder_addendum):

    import mysql.connector
    import loompy
    
    from selenium.common.exceptions import NoSuchElementException, TimeoutException
    import json

    import time

    with open('/etc/dt-config.json') as config_file:
        config = json.load(config_file)

    import sys
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
        
        os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
        
        sources = context["Sources_global"]["Sources"]
        orthology_sources = ["MozAtlas", "MozTubules", "AegyptiAtlas", "SilkDB"]
        orthology_species = {
            'Hsap': 'Homo sapiens',
            'Agam':'Anopheles gambiae',
            'Aaeg':'Aedes aegypti', 
            'Bmor':'Bombyx mori'
        }

        #Generic settings for connecting to MySQL
        #This may be removed in the future if SQL migration is necessary
        sql_settings = {
            'host':config['SQL_HOST'],
            'user':config['SQL_USER'],
            'pass':config['SQL_PASS']
        }

        extra_threshold_sources = []
        
        #Handles User-Uploaded datasets, passes needed variables to the Agglomerator
        if "UsrUpload_SETTINGS" in context:
            print("Processing user-uploaded data!")
            usr_metricweights = f'GLOBAL:{return_global_weight(context, "UserData")}'
            print("Weights loaded")
            usr_metricweights, usr_location = UserUpload(context['UsrUpload_SETTINGS'], context['UploadedData'], usr_metricweights, outputfolder, config)
        else:
            usr_location = ''
            usr_metricweights = 'GLOBAL:0;Enrichment:0;Specificity:0'

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

        try:
            print('Initialising the Tally')
            StartTallying(context, 
                          fa1_location, 
                          fa1_extra_location, 
                          fa1weightstring, 
                          fa2_location, 
                          fa2_extra_location, 
                          fa2weightstring, 
                          outputfolder,
                          usr_location,
                          usr_metricweights)
            
        except Exception as e:
            print(e)

        if "FlyCellAtlas" in sources:
            print('Processing FlyCellAtlas!')
            if 'References' in context:
                ref_genes = context['References']
            else:
                ref_genes = ''

            fca_metricweights = f'GLOBAL:{return_global_weight(context, "FlyCellAtlas")}'
            
            loom_inc = 0
            finished = 0

            import time

            while loom_inc < 10 and finished == 0:
                loomfile = f'AssociatedFiles/FCA_WholeFly_Loom_{loom_inc}.loom'

                try:
                    fca_metricweights = FCA(context['FlyCellAtlas'], ref_genes, fca_metricweights, outputfolder,loomfile)
                    finished = 1

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

        if "FlyBase" in sources:
            print('Processing FlyBase!')
            
            os.mkdir(f'{outputfolder}/FlyBase_Analysis')
            os.mkdir(f'{outputfolder}/FlyBase_Analysis/Downloads')
            os.mkdir(f'{outputfolder}/FlyBase_Analysis/Downloaded_lists')
            
            done = 0
            
            while done == 0:
                try:
                    fb_metricweights = f'GLOBAL:{return_global_weight(context, "FlyBase")}'
            
                    fb_metricweights = FB(context['FlyBase'], fb_metricweights, outputfolder)
                    done = 1
                    
                except NoSuchElementException as e1:
                    print(e1)
                    print('\nLIKELY CAUSE: FLYBASE LAYOUT CHANGE\n')
                    fb_metricweights = 'GLOBAL:0;Association:0;Phenotype:0'
                    sources.remove('FlyBase')
                    done=1
                except TimeoutException as e2:
                    print(e2)
                    print('\nLIKELY CAUSE: CONNECTION ISSUES ON SERVER SIDE\n')
                    fb_metricweights = 'GLOBAL:0;Association:0;Phenotype:0'
                    sources.remove('FlyBase')
                    done=1
                except ValueError as e3:
                    print(e3)
                    print('\nLIKELY CAUSE: REPEATRED FAILURE \n')
                    fb_metricweights = 'GLOBAL:0;Association:0;Phenotype:0'
                    sources.remove('FlyBase')
                    done=1

            print('\n FLYBASE DONE - %s seconds' % (time.time() - last_time))
            last_time = time.time()

        else:
            fb_metricweights = 'GLOBAL:0;Association:0;Phenotype:0'

        species = []

        if 'Orthology' in context['Scoring']['Metrics']:
            
            speciesweightstring = f'GLOBAL:{context["Scoring"]["Weights"]["Orthology"]}'

            os.mkdir(f'{outputfolder}/Orthology')
            
            tmpspecies = []
            
            for ortho_species in orthology_species:

                species_weight = context["Scoring"]["Weights"][ortho_species]
                speciesweightstring = speciesweightstring + f';{ortho_species}:{species_weight}'
                
                if context["Scoring"]["Weights"][ortho_species] != 0:
                    tmpspecies.append(orthology_species[ortho_species])
            
            print(tmpspecies)
            
            if "DIOPT" in sources:
                print('Processing DIOPT!')
                species.append('Homo sapiens')

                HSAP(context['DIOPT'], outputfolder)

                print('\n DIOPT DONE - %s seconds' % (time.time() - last_time))
                last_time = time.time()

            if any(item in sources for item in orthology_sources):

                if 'Homo sapiens' in tmpspecies:
                    tmpspecies.remove('Homo sapiens')

                print('Finding Orthologs!')
                
                ORTHODB(tmpspecies, outputfolder)

                print('\n ORTHOLOG ID DONE - %s seconds' % (time.time() - last_time))
                last_time = time.time()

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

            if "AegyptiAtlas" in sources:
                print('Processing AegyptiAtlas!')
                species.append('Aedes aegypti')
                AAEG(context['AegyptiAtlas'], outputfolder)

                print('\n AAEG DONE - %s seconds' % (time.time() - last_time))
                last_time = time.time()

            if "SilkDB" in sources:
                print('Processing SilkDB!')
                species.append('Bombyx mori')
                BMOR(context['SilkDB'], outputfolder)

                print('\n BMORI DONE - %s seconds' % (time.time() - last_time))
                last_time = time.time()
        else:
            speciesweightstring = 'GLOBAL:0;Hsap:0;Agam:0;Aaeg:0;Bmor:0'

        global_score_weights = gather_global_weights(context)
        
        print('Finishing the Analysis!')
        print(fb_metricweights)

        try:
            CompleteAnalysis(
                global_score_weights,
                usr_metricweights, 
                fa1_metricweights, 
                fa2_metricweights, 
                fca_metricweights, 
                fb_metricweights, 
                speciesweightstring, 
                species, 
                outputfolder)
        except Exception as e:
            print(e)

        print('\n FINAL TALLYING DONE - %s seconds' % (time.time() - last_time))
        print('\n\n DIGITtally DONE - %s seconds TOTAL' % (time.time() - start_time))

