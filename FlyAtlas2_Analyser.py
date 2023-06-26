#Code written by Dr. Andrew Gillen (Dow-Davies lab, University of Glasgow)
#Originally written October 2022, last updated 26/06/2023

import csv
import argparse
import os

import mysql.connector
from mysql.connector import Error

#Part of the pipeline for DIGITtally analysis - see www.digittally.org
#Specifically, this program analyses the microarray data present in FlyAtlas2 (https://motif.mvls.gla.ac.uk/FlyAtlas2/)

#Connects to the specified MySQL database
def create_server_connection(host_name, user_name, user_password, db_name):
    connection = None
    try:
        connection = mysql.connector.connect(
            host=host_name,
            user=user_name,
            passwd=user_password,
            database=db_name
        )
        print("MySQL Database connection successful")
    except Error as err:
        print(f"Error: '{err}'")

    return connection

#Submits a MySQL query
def read_query(connection, query):
    cursor = connection.cursor()
    result = None
    try:
        cursor.execute(query)
        result = cursor.fetchall()
        return result
    except Error as err:
        print(f"Error: '{err}'")

#Builds dictionaries associating codes from SQL database to Fly Type and Tissue.
def populate_tissues(connection):
    tissuequery = """
    SELECT *
    FROM Tissue;
    """
    tissuedata = read_query(connection, tissuequery)

    flycat = {}
    tissuecat = {}

    for line in tissuedata:
        if line[1] == 'Larval':
            flycat[line[0]] = line[1]
        else:
            flycat[line[0]] = line[2]

        tissuecat[line[0]] = line[3]
    
    return(flycat, tissuecat)

#Defines background levels (based on FlyAtlas2 specifications) and ID formats, allowing analysis at both the compiled gene and transcript level, for protein coding and miRNA genes.
def get_method(mode, supplied_fpkm_thresh):

    if mode <= 2:
        quanttype = 'FPKM'
        threshold = supplied_fpkm_thresh
    elif mode > 2:
        threshold = 100.0
        quanttype = 'RPM'
    
    if mode % 2 == 1:
        target = 'Transcript'
        idformat = 'FBtr'
    else:
        target = 'Gene'
        idformat = 'FBgn'
    
    return(quanttype, threshold, target, idformat)

#Defines "whole fly" reads from each fly type, used to calculate enrichment. This is adjusted by the weights calculated by TGF_Weightfinder
def populate_wholefly(connection, flytypes, methodcode, weights, unweighted, supplied_fpkm_thresh):
    wholeflydict = {}
    passlist = []

    for type in flytypes:
        wholeflydict[type] = {}
    
    quanttype, threshold, target, idformat = get_method(methodcode, supplied_fpkm_thresh)

    #MySQL Query for MULTIPLE tables in the FlyAtlas2 database-(FPKM/RPM)(Gene/Transcript). Genes with a status other than OK (from RNAseq quantification output) are not included. 
    #Two queries created as RPM tables do not have a status column
    if methodcode <=2:
        wholefly_query = f"""
        SELECT {idformat}, TissueID, {quanttype}, Status
        FROM {target}{quanttype}
        WHERE TissueID IN (100, 200, 300)
        """

    else:
        wholefly_query = f"""
        SELECT {idformat}, TissueID, {quanttype}
        FROM {target}{quanttype}
        WHERE TissueID IN (100, 200, 300)
        """
    
    print('Whole Fly query: ', wholefly_query)
    wholeflyresults = read_query(connection, wholefly_query)

    print('query successful')

    #Adjusts FPKM/RPM readings based on Median of Ratio weights. Only weights when "unweighted" arg is not altered
    if unweighted != 0:
        with open(f'{weights}/{target}{quanttype}.csv') as weightfile:
            weightsIN = csv.reader(weightfile)

            for line in weightsIN:
                if line[0] == 'Whole body':
                    wbm = float(line[1])
                    wbf = float(line[2])
                    wbl = float(line[3])

    for line in wholeflyresults:
        readid = line[0]
        tissueid = line[1] 
        read_value = float(line[2])

        try:
            status = line[3]
        except:
            status = 'OK'

        if status == 'OK' :
            #Given that the background level of detection should be the absolute minimum value used to calculate enrichment, to prevent skewing results due to abnormally low "Whole Fly" readings,
            #Similar blocks like this are used throughout to set the minimal level of "Whole Fly" expression to the appropriate threshold (2.0 FPKM, 100 RPM)
            if unweighted == 0:
                read_value = max([read_value, threshold])
            
            #TissueID 100 in the FlyAtlas2 database is Male Whole Fly
            if tissueid == 100:
                if unweighted != 0:
                    read_value = max([read_value/wbm, threshold])

                wholeflydict['Male'][readid] = read_value

            #TissueID 200 in the FlyAtlas2 database is Female Whole Fly
            elif tissueid == 200:
                if unweighted != 0:
                    read_value = max([read_value/wbf, threshold])

                wholeflydict['Female'][readid] = read_value
            
            #TissueID 300 in the FlyAtlas2 database is Larval Whole Fly
            elif tissueid == 300:
                if unweighted != 0:
                    read_value = max([read_value/wbl, 2.0])

                wholeflydict['Larval'][readid] = read_value

        #This block is necessary to pass over genes which have an exit condition other than OK in later analyses
        else:
            if tissueid == 100:
                failtype = 'Male'
            elif tissueid == 200:
                failtype = 'Female'
            elif tissueid == 300:
                failtype = 'Larval'
            
            passlist.append((readid, failtype))

    return(wholeflydict, passlist)

#Gathers Abundance (Raw FPKM/RPM) and Enrichment (FPKM/Whole Fly) readings for each tissue for each fly type. 
#These are adjusted by the weights calculated by TGF_Weightfinder
def get_expression(connection, flycatdict, tissuecatdict, wholeflydict, flytypes, methodcode, weights, unweighted, passlist, supplied_fpkm_thresh):

    quanttype, threshold, target, idformat = get_method(methodcode, supplied_fpkm_thresh)

    #MySQL Query for MULTIPLE tables in the FlyAtlas2 database-(FPKM/RPM)(Gene/Transcript)
    query = f"""
    SELECT {idformat}, TissueID, {quanttype}
    FROM {target}{quanttype}
    """

    expressionresults = read_query(connection, query)
    weightbytisstype = {}

    resultsdict = {}
    abundancedict = {}

    if unweighted != 0:
        with open(f'{weights}/{target}{quanttype}.csv') as weightfile:
            weightsIN = csv.reader(weightfile)

            for line in weightsIN:
                tiss = line[0]
                if tiss == 'Brain' or tiss == 'CNS':
                    tiss = 'Brain/CNS'

                try: 
                    tiss_m = float(line[1])
                except:
                    tiss_m = 1
                try:
                    tiss_f = float(line[2])
                except:
                    tiss_f = 1
                try:
                    tiss_l = float(line[3])
                except:
                    tiss_l = 1

                weightbytisstype[f'Male_{tiss}'] = tiss_m
                weightbytisstype[f'Female_{tiss}'] = tiss_f
                weightbytisstype[f'Larval_{tiss}'] = tiss_l

    for line in expressionresults:
        readid = line[0]

        if readid not in resultsdict:
            resultsdict[readid] = {}
            abundancedict[readid] = {}

            for type in flytypes:
                resultsdict[readid][type] = {}
                abundancedict[readid][type] = {}

        tissueid = line[1]
        read_value = line[2]
        flytype = flycatdict[tissueid]
        tissuetype = tissuecatdict[tissueid]

        if tissuetype == 'Brain' or tissuetype == 'CNS':
            tissuetype = 'Brain/CNS'
        
        #excludes reads below the set background threshold 
        if read_value < threshold:
            read = float('nan')

        else:
            #Only weights when "unweighted" arg is not altered
            if unweighted != 0:
                read_value = (read_value/weightbytisstype[f'{flytype}_{tissuetype}'])

            if ((readid,flytype) in passlist): 
                read = float('nan')
            elif read_value < threshold:
                read = 0
            else:
                read = read_value/wholeflydict[flytype][readid]

        resultsdict[readid][flytype][tissuetype] = read
        abundancedict[readid][flytype][tissuetype] = read_value

    return resultsdict, abundancedict

#Gathers tissues of interest from designated text files
def build_tiss_lists(tissuesfile):
    fulllist = []

    with open(tissuesfile) as INtissues:
        for line in INtissues:
            tissue = line.strip('\n')
            fulllist.append(tissue)

    return fulllist

#Generates output files, both for further pipeline processing and for general analysis
def build_output(methodcode, enrichmentdict, abundancedict, flytypes, enrstringency, abnstringency, toi, outfolder, permissive, supplied_fpkm_thresh):


    quanttype, threshold, target, idformat = get_method(methodcode, supplied_fpkm_thresh)
    
    #Sets up dictionaries necessary for building outputs
    toidict = {}
    allbutdict = {}
    toi_enrich_dict = {}
    toi_abun_dict = {}
    enriched_genes = {}
    abundant_genes = {}
    which_tissues_enr = {}
    which_tissues_abn = {}

    for tissue in toi:

        toidict[tissue] = {}
        allbutdict[tissue] = {}

        for type in flytypes:
            toidict[tissue][type] = {}
            allbutdict[tissue][type] = {}
            toi_enrich_dict[type] = {}
            toi_abun_dict[type] = {}
            enriched_genes[type] = []
            abundant_genes[type] = []
            which_tissues_enr[type] = {}
            which_tissues_abn[type] = {}

    for type in flytypes:
        referenceheader = [f'{quanttype}', ]
        headerprinted = 0

        #This block gathers info necessary for the HitScore summary, whilst also generating four .csvs, 
        #Two "tissuefiles": one with info on genes/transcripts ENRICHED in at least one tissue, the other with info on increased ABUNDANCE in tissues of interest
        #Two "referencefiles" for ease of checking data on interesting genes - ie those enriched/increased abundance in at least one tissue of interest
        with open(f'{outfolder}/{type}_{quanttype}{target}_EnrichedByTissue.csv', 'w+') as tissuefileenr, open(f'{outfolder}/{type}_{quanttype}{target}_AbundantByTissue.csv', 'w+') as tissuefileabn, open(f'{outfolder}/{type}_{quanttype}{target}_AbundanceReference.csv', 'w+') as referencefileabn, open(f'{outfolder}/{type}_{quanttype}{target}_EnrichmentReference.csv', 'w+') as referencefileenr:

            OUTtissue_enr = csv.writer(tissuefileenr)
            OUTtissue_abn = csv.writer(tissuefileabn)
            OUTref_enr = csv.writer(referencefileenr)
            OUTref_abn = csv.writer(referencefileabn)

            for product in enrichmentdict:

                #referencefile line only to be printed if the gene looks interesting
                refline_enr = [product]
                refline_abn = [product]
                printref = 0

                toi_enrichments = []
                toi_abundances = []
                nontarget_abundances = []

                which_tissues_enr[type][product] = []
                which_tissues_abn[type][product] = []

                #All other tissues must be checked against every tissue of interest - the "allbutdict" allows this
                for tissueofinterest in toi:
                    allbutdict[tissueofinterest][type][product] = []
                
                toi_enrich_dict[type][product] = 0

                #gathers information on each gene for each tissue
                for tissue in enrichmentdict[product][type]:
                    if tissue not in referenceheader:
                        referenceheader.append(tissue)

                    enrichment = enrichmentdict[product][type][tissue]
                    abundance = abundancedict[product][type][tissue]

                    refline_enr.append(enrichment)
                    refline_abn.append(abundance)
                    
                    #For the tissues of interest, enrichment > threshold is recorded
                    if tissue in toi:
                        
                        if enrichment >= enrstringency:
                            toi_enrich_dict[type][product] += 1
                            which_tissues_enr[type][product].append(tissue)

                        toi_enrichments.append(enrichment)
                        toi_abundances.append(abundance)

                        toidict[tissue][type][product] = enrichment

                        for othertissue in toi:
                            if tissue != othertissue:
                                allbutdict[othertissue][type][product].append(enrichment)

                    #Other tissues are just used to define non target enrichment levels
                    else:
                        if tissue != 'Whole body' and tissue not in permissive:
                            nontarget_abundances.append(abundance)

                            for othertissue in toi:
                                if tissue != othertissue:
                                    allbutdict[othertissue][type][product].append(enrichment)

                tissuehead = ['ID', 'Tissues', 'No. Enriched Tissues of Interest']

                toi_abun_dict[type][product] = 0
                
                #If a gene is enriched in ALL tisses of interest, it will be printed to the final Enrichment Genelist
                if toi_enrich_dict[type][product] == len(toi):
                    enriched_genes[type].append(product)
                
                #Checks which tissues of interest express a gene over (abundance stringency * highest non-target tissue)
                for tissue in toi:
                    try:
                        if round(abundancedict[product][type][tissue]) > round(abnstringency * max(nontarget_abundances)):
                            toi_abun_dict[type][product] += 1
                            which_tissues_abn[type][product].append(tissue)
                    except:
                        pass
                
                #If a gene is increased in abundance in ALL tisses of interest, it will be printed to the final Abundance Genelist
                if toi_abun_dict[type][product] == len(toi):
                    abundant_genes[type].append(product)
                                
                #prints headers for each output file
                if headerprinted == 0:
                    OUTtissue_enr.writerow(tissuehead)
                    OUTtissue_abn.writerow(tissuehead)
                    OUTref_enr.writerow(referenceheader)
                    OUTref_abn.writerow(referenceheader)
                    headerprinted += 1

                #Checks whether gene is enriched in any tissue of interest
                if toi_enrich_dict[type][product] > 0:
                    compartoutrow = [product, which_tissues_enr[type][product], len(which_tissues_enr[type][product])]
                    OUTtissue_enr.writerow(compartoutrow)
                    printref = 1
                
                #Checks whether gene is at increased abundance in any tissue of interest
                if toi_abun_dict[type][product] > 0:
                    compartoutrow = [product, which_tissues_abn[type][product], len(which_tissues_abn[type][product])]
                    OUTtissue_abn.writerow(compartoutrow)
                    printref = 1
                
                #Prints the reference info for interesting genes
                if printref == 1:
                    OUTref_enr.writerow(refline_enr)
                    OUTref_abn.writerow(refline_abn)

        #This block creates the final Enriched genelists - 
        #First individual files for each fly type are generated
        for indivtype in flytypes:
            with open(f'{outfolder}/{quanttype}{target}_ENRICHEDin{indivtype.upper()}.txt', 'w+') as indivtypefile:
                for enrgene in enriched_genes[indivtype]:
                    indivtypefile.write(f'{enrgene}\n')

        #Next, genes which are enriched in all adult flies (ie male + female) are printed to one list, and genes enriched in ALL fly types to another.
        with open(f'{outfolder}/{quanttype}{target}_ENRICHEDinADULTS.txt', 'w+') as malefemalefile, open(f'{outfolder}/{quanttype}{target}_ENRICHEDinALL.txt', 'w+') as larvalfile, open(f'{outfolder}/{quanttype}{target}_ENRICHEDinALL.txt', 'w+') as alltypesfile:
            male_female = set(enriched_genes['Male']) & set(enriched_genes['Female'])

            for item in male_female:
                malefemalefile.write(f'{item}\n')

            all_types = set(male_female) & set(enriched_genes['Larval'])

            for item in all_types:
                alltypesfile.write(f'{item}\n')

        #This block creates the final Abundance genelists - 
        #First individual files for each fly type are generated
        for indivtype in flytypes:
            with open(f'{outfolder}/{quanttype}{target}_ABUNDANTin{indivtype.upper()}.txt', 'w+') as indivtypefile:
                for abngene in abundant_genes[indivtype]:
                    indivtypefile.write(f'{abngene}\n')

        #Next, genes which are increased in abundance in all adult flies (ie male + female) are printed to one list, and genes with increased abundance in ALL fly types to another.
        with open(f'{outfolder}/{quanttype}{target}_ABUNDANTinADULTS.txt', 'w+') as malefemalefile_abn, open(f'{outfolder}/{quanttype}{target}_ABUNDANTinALL.txt', 'w+') as alltypesfile_abn:
            male_female = set(abundant_genes['Male']) & set(abundant_genes['Female'])

            for item in male_female:
                malefemalefile_abn.write(f'{item}\n')

            all_types = set(male_female) & set(abundant_genes['Larval'])

            for item in all_types:
                alltypesfile_abn.write(f'{item}\n')

#Converts the float values for thresholds to strings, replacing "."s with "-"s for filenames
def make_string(floating_point, identifier):
    string_vers = str(floating_point)
    string_vers = string_vers.replace('.', '-')

    final_string = identifier + string_vers

    return final_string

def analyse_FA2(database, host, pw, user, enrstringency, abnstringency, weights, weightingmode, targettissues, permissive, fpkmthresh, outfolder):

    flytypes = ['Male', 'Female', 'Larval']

    connection = create_server_connection(host, user, pw, database)

    flycat, tissuecat = populate_tissues(connection)
    print('Data populated')

    toi = targettissues

    for i in range(1,5):
        wholefly, passlist = populate_wholefly(connection, flytypes, i, weights, weightingmode, fpkmthresh)
        print('wholefly setup ok')

        enrichmentresults, abundanceresults = get_expression(connection, flycat, tissuecat,  wholefly, flytypes, i, weights, weightingmode, passlist, fpkmthresh)
        print('expression data gathered ok')

        build_output(i, enrichmentresults, abundanceresults, flytypes, enrstringency, abnstringency, toi, outfolder, permissive, fpkmthresh)
        print('output generated')
