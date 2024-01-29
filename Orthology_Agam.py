import csv

from collections import defaultdict
import re

import mysql.connector
from mysql.connector import Error

#A simple function to break up the weight string supplied in input arguments into an accessible dictionary
def gather_weights(targetweights):
    weight_dict = {}

    weights = targetweights.split(';')

    #parses the supplied weight into a dictionary
    for item in weights:
        subitem = item.split(':')
        indivweight = float(subitem[1])
        weight_dict[subitem[0]] = indivweight

    return weight_dict

def define_sexes(target_dict):
    target_sexes = set()

    if target_dict['MALE'] > 0:
        target_sexes.add('Male')
    
    if target_dict['FEMALE'] > 0:
        target_sexes.add('Female')
    
    if target_dict['ADULT'] > 0:
        target_sexes.add('Male')
        target_sexes.add('Female')
    
    return target_sexes

#Gathers tissues of interest from designated text files
def build_tiss_lists(tissuesfile):
    fulllist = set()

    with open(tissuesfile) as INtissues:
        for line in INtissues:
            tissue = line.strip('\n')
            fulllist.add(tissue)

    return fulllist

#Creates the default tally for MozAtlas
def mozatlas_default():
    return [0, 0, 0]

#Creates the default tally for MozTubule
def moztubule_default():
    return [0, 0]

#Populates the A. gam orthologs for each gene of interest, using the designated ortholog file
def pop_orthos(input):
    ortholist = set()
    problem_orthos = {}

    dmeltoorthodict = defaultdict(list)
    orthotodmeldict = {}

    has_an_ortho = []
    all_genes = []

    with open(input) as foundorthosfile:
        INfos = csv.reader(foundorthosfile)

        next(INfos)

        for line in INfos:
          
            fbgene = line[0]
            all_genes.append(fbgene)

            if line[1] != 'No orthologs found':

                #if more than one ortholog is present, these are split and held in a list, to be processed seperately
                if ';' in line[1]:
                    indiv_orthos = line[1].split('; ')
                else:
                    indiv_orthos = [line[1]]

                has_an_ortho.append(fbgene)

                for gene in indiv_orthos:
                    gene = gene.strip(' ')
                    #we set up a dictionary of lists which holds all the orthologs associated with each FbID
                    dmeltoorthodict[fbgene].append(gene)

                    #The reverse is slightly harder - we want to create a dictionary matching each ortholog to a specific FbID but some orthologs can be associated with multiple
                    if gene not in ortholist:
                        orthotodmeldict[gene] = fbgene
                        ortholist.add(gene)

                    #A.gam genes which are already associated with a FbID are added to the "problem ortho" dictionary
                    else:
                        
                        #This ensures correct layout of the dict - {A.gam gene: [list of flybase IDs]}
                        try:
                            problem_orthos[gene].append(fbgene)
                        except:
                            problem_orthos[gene] = [orthotodmeldict[gene], fbgene]
    print('finished function')
    return ortholist, orthotodmeldict, dmeltoorthodict, problem_orthos, has_an_ortho, all_genes

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

#Uses the MozAtlas "gene_table" table to build a list of probe IDs associated with each ortholog of interest
def get_probe_ids(connection, orthologs):

    probeids = defaultdict(set)

    for indivortholog in orthologs:
        query = f"""SELECT gene_id, probe_id
        FROM `mozatlas_gene_table`
        WHERE `gene_id` = '{indivortholog}'
        """

        probeidqres= read_query(connection, query)

        if probeidqres != []:
            for item in probeidqres:
                if item[1] != '':
                    probeids[indivortholog].add(item[1])    

    print('probe ids OK')             
    return probeids

#Gets the tissue and sex specific expression of a given probe
def get_tissue_expression(probe, connection):
    tissquery = """SELECT *
                FROM `mozatlas_expression`
                WHERE `probe_id` = '{_probe}'
                """.format(_probe = probe)

    tq = read_query(connection, tissquery)


    if tq == []:
        return "no data"
    
    return tq[0][1:]

#checks whether expression meets the presence threshold - if so, values are calculated normally
#readval = average expression, enrval = avg expression/whole body expression
def check_presence(presence_threshold, tissuedata, wholebody):

    if tissuedata['Present'] < presence_threshold:
        readval = 0
        enrval = 0

    else:
        readval = tissuedata['Mean']
        enrval = (readval / wholebody)
    
    return readval, enrval

#Sets up the following system for marking off individual tissues of interest in which orthologs of a given gene have been found:
    #{Gene:
    #   {Sex:       
    #       {Category:
    #           [Individual, Tissues]
    #       }
    #   }
    # }
#By removing tissues one at a time when any ortholog of a given Drosophila gene fulfils a criteria, that "category" dictionary will be empty if necessary conditions are fulfilled in all tissues of interest
#Thus, we can track, across various orthologs, many of which have multiple Microarray probes associated, which Drosophila genes have orthologs Enriched/Abundant/Expressed in a set of tissues of interest
def initialise_scoring_system(has_an_ortho, targettiss, sexes):
    uncompleted_cats = {}

    for gene in has_an_ortho:
        uncompleted_cats[gene] = {}

        for sex in sexes:
            uncompleted_cats[gene][sex] = {}

            for category in ['Enrichment', 'Specificity', 'Expression']:
                uncompleted_cats[gene][sex][category] = [] + list(targettiss)
    
    return uncompleted_cats

#Creates a dictionary containing details of tissues to be removed from the uncompleted_cats scoring system
#{sex:
#   {category:
#       [Individual, Tissues]
#   }
#}
def removal_by_ortholog(sexes):
    to_remove = {}

    for sex in sexes:
        to_remove[sex] = {}

        for category in ['Enrichment', 'Specificity', 'Expression']:
            to_remove[sex][category] = []
    
    return to_remove

#Uses the "to_remove" dictionary to remove individual tissues for which conditions have been fulfilled from the "Uncompleted_cats" dictionary on a BY GENE basis
#Notably, this can be called on genes which have already been processed previously, removing tissues in which conditions may not have been met on the previous pass
#THIS IS A VITAL FEATURE
def effect_removal(to_remove, uncompleted_cats_section):
    newcats = uncompleted_cats_section
    
    for sex in to_remove:

        for category in to_remove[sex]:

            for tissue in to_remove[sex][category]:

                if tissue in newcats[sex][category]:

                    newcats[sex][category].remove(tissue)
    
    return newcats

#This function was designed to speed up MozAtlas Analysis - by gathering all tissues at once we reduce the number of necessary SQL calls.
def find_mozatlas_data_cols(connection):
    print('Starting column processing')
    def dictdict():
        return defaultdict(dict)

    indices = defaultdict(dictdict)

    #This function was designed to speed up MozAtlas Analysis - by gathering all tissues at once we reduce the number of necessary SQL calls.
    colquery = """SELECT COLUMN_NAME, ORDINAL_POSITION
            FROM INFORMATION_SCHEMA.COLUMNS
            WHERE TABLE_NAME = N'mozatlas_expression'
            """

    cq = read_query(connection, colquery)

    for data_segment in cq:

        try:
            parsed_word = re.findall('[a-zA-Z][^A-Z]*', data_segment[0])
            
            if len(parsed_word) > 3:
                parsed_word[2] = "".join(parsed_word[2:])

            indices[parsed_word[0]][parsed_word[1]][parsed_word[2]] = (data_segment[1]-2)

        except:
            pass

    print('Columns OK')
    return indices

#This code handles actually analysing MozAtlas data, using the uncompleted_cats system to track the genes for which specific tissues have fulfilled specific criteria.
def find_mozatlas_expression(connection, has_an_ortho, target_sex, targettiss, probeids, orthotodmel, enrichment_threshold, specificity_threshold, presence_threshold, problemorthos):
    
    #This block sets up static lists, defining the tissues available in MozAtlas.
    #The first list holds tissues with data for both sexes
    tissues = [('Carcass', 'Carcass'), ('Head','Head'), ('Malpighian', 'Malpighian Tubule'), ('Midgut', 'Midgut'), ('Salivary', 'Salivary Gland')]

    #This sexspecific_tiss dictionary holds a male and female list, with tissues which are sex-specific and thus only have one set of data associated
    sexspecific_tiss = {}
    sexspecific_tiss['male'] = [('Acps', 'Male Accessory Glands'), ('Testis', 'Testis')]
    sexspecific_tiss['female'] = [('Ovary', 'Ovary')]

    #This short list just holds the two categories of data we need from MozAtlas - Gene mean expression and Gene present counts
    datacat = ['Mean', 'Present']

    #Processed the user-input desired sexes into usable format
    sexes = [x.lower() for x in target_sex]

    #Sets up uncompleted_cats
    uncompleted_cats = initialise_scoring_system(has_an_ortho, targettiss, sexes)

    ma_indices = find_mozatlas_data_cols(connection)

    #Each A.gam gene is processed seperately
    for mozgene in probeids:

        to_remove = {}
        
        #Creates a list in "to_remove" for each sex of interest
        for sex in sexes:
            to_remove[sex] = defaultdict(list)

        #Gets the flybase gene associated with the A.gam gene being processed.
        fbgene = orthotodmel[mozgene]

        #Each individual microarray probe is handled seperately
        for probe in probeids[mozgene]:
            probe_expression = {}

            mozatlas_data = get_tissue_expression(probe, connection)


            if mozatlas_data != 'no data':

                for sex in sexes:
                    sexnontargets = []

                    #The whole body reading is processed first to calculate enrichments
                    try:
                        wholebody = float(mozatlas_data[ma_indices[sex]['Body']['Mean']])

                    except Error as e:
                        print(e)
                        print(probe, 'STILL CAUSING CRASHES')
                        exit()

                    #Tissues present in both sexes are processed first
                    for giventissue in tissues:
                        tissuedata = {}

                        #Mean expression and present call values are gathered
                        for datatype in datacat:

                            tissuedata[datatype] = float(mozatlas_data[ma_indices[sex][giventissue[0]][datatype]])

                        
                        #The present call values are used to adjust enrichment/abundance values
                        readval, enrval = check_presence(presence_threshold, tissuedata, wholebody)
                        
                        if giventissue[1] in targettiss:
                            probe_expression[giventissue[1]] = readval

                            #Tissues of interest enriched over the threshold are marked in the Enrichment "to_remove" list
                            if enrval > enrichment_threshold:
                                to_remove[sex]['Enrichment'].append(giventissue[1])

                            #Tissues of interest showing any expression are marked in the Abundance "to_remove" list
                            if readval > 0:
                                to_remove[sex]['Expression'].append(giventissue[1])

                        #Non-target tissue read values are added to the sex-specific "sexnontargets" list
                        else:
                            sexnontargets.append(readval)
                    
                    #The sex-specific tissues for the sex currently being processed are handled next
                    for giventissue in sexspecific_tiss[sex]:
                        
                        #Mean expression and present call values are gathered
                        for datatype in datacat:

                            tissuedata[datatype] = float(mozatlas_data[ma_indices[sex][giventissue[0]][datatype]])

                        #The present call values are used to adjust enrichment/abundance values
                        readval, enrval = check_presence(presence_threshold, tissuedata, wholebody)

                        if giventissue[1] in targettiss:
                            probe_expression[giventissue[1]] = readval
                            
                            #Tissues of interest enriched over the threshold are marked in the Enrichment "to_remove" list
                            if enrval > enrichment_threshold :
                                to_remove[sex]['Enrichment'].append(giventissue[1])

                            #Tissues of interest showing any expression are marked in the Abundance "to_remove" list
                            if readval > 0:
                                to_remove[sex]['Expression'].append(giventissue[1])
                        
                        #Non-target tissue read values are added to the sex-specific "sexnontargets" list
                        else:
                            sexnontargets.append(readval)

                    #Now that all non-target tissue read values have been added to sexnontargets, we can check, for each tissue of interest, whether expression is higher than ALL non-target tissues.
                    for tissue in targettiss:
                        if probe_expression[tissue] > (max(sexnontargets) * specificity_threshold):
                            to_remove[sex]['Specificity'].append(tissue)

        #If the A.gam gene is an ortholog of multiple D.mel genes, the fulfillment of conditions from tissues of interest is recorded for ALL the appropriate D.mel genes
        if mozgene in problemorthos:
            for indiv_fbid in problemorthos[mozgene]:
                uncompleted_cats[indiv_fbid] = effect_removal(to_remove, uncompleted_cats[indiv_fbid])
        
        #If the A.gam gene is an ortholog of only one D.mel gene, the fulfilment of conditions is recorded only for this gene
        else:
            uncompleted_cats[fbgene] = effect_removal(to_remove, uncompleted_cats[fbgene])
    
    print('MozAtlas Done')
    return uncompleted_cats

#The uncompleted_cats dictionary is resolved by counting empty lists
#If the number of empty lists for each gene, for each category, is equal to the number of sexes of interest, that category scores a "1" in the MozAtlas tally
#Otherwise, the category scores 0
def score_by_cat(uncompleted_cats, weight_dict):

    tally_compilation = {}
    
    pos_by_cat = {
        'Enrichment': 0, 'Specificity': 1, 'Expression': 2
    }

    for flytype in ['MALE', 'FEMALE', 'ADULT']:
        if weight_dict[flytype] != 0:
            tally_compilation[flytype] = defaultdict(mozatlas_default)
        
    for gene in uncompleted_cats:
        
        gene_cats = {'Enrichment': 0, 'Specificity': 0, 'Expression': 0}
        num_sexes = 0

        for sex in uncompleted_cats[gene]:
            num_sexes += 1

            for category in uncompleted_cats[gene][sex]:

                if uncompleted_cats[gene][sex][category] == []:
                    gene_cats[category] += 1
                    
                    sex_for_comp = sex.upper()
                    tally_compilation[sex_for_comp][gene][pos_by_cat[category]] = 1

        
        #The ADULT flytype only scores if both male and female sexes score
        if 'ADULT' in tally_compilation:

            for cat in pos_by_cat:

                if gene_cats[cat] == num_sexes:
                    tally_compilation['ADULT'][gene][pos_by_cat[cat]] = 1

    return tally_compilation

#Moztubule uses a seperate microarray schema, so probe ids must be gathered seperately for this system.
def get_moztubule_probeids(orthologs, probefile):

    probeids = defaultdict(set)

    genetoprobes = defaultdict(set)

    with open(probefile) as mtannofile:
        INmtannos = csv.reader(mtannofile)

        for i in range(23):
            next(INmtannos)

        #The positioning of gene names is inconsistent - this block gathers the name wherever it may be found.
        for line in INmtannos:

            gene = line[17]

            if gene == '---':
                gene = line[13]
                gene = gene[:-3]

            if gene[0:4] != 'AGAP':
                gene = float('nan')

            probe = line[0]

            genetoprobes[gene].add(probe)

    #The probe IDs corresponding to A.gam orthologs to D.mel genes of interest are stored in probeids.
    for indivortho in orthologs:
        try:
            for probe in genetoprobes[indivortho]:
                probeids[indivortho].add(probe)

        except:
            pass

    return probeids

#This function parses the expression data found in MozTubules, creating a dictionary (mtdata) which holds expression data for each gene (key = Agam gene symbol)
def extract_mt_expression(expressionfile):
    mtdata = {}

    with open(expressionfile) as mtdatafile:

        linecount = 0

        for line in mtdatafile:
            if linecount == 0:
                linecount += 1
            
            else:
                parsedline = line.split('\t')
                readinfo = []

                for i in range(1, 7):
                    readinfo.append(parsedline[i])

                mtdata[parsedline[0]] = readinfo
    
    return mtdata

#This function gathers the expression of each Agam gene of interest from the MozTubules dataset for all flytypes of interest
#MozTubules has a seperate "Adult" score - hence why adult may find a hit here but either male of female may fail the enrichment threshold.
def find_moztubules_expression(moztubprobeids, orthotodmel, threshold, expressionfile, problemorthos, weight_dict, has_an_ortho):
    
    moztubule_tally = {}
    true_moztubule_tally = {}

    default_gene_score = []
    true_default = []

    pos_in_weights = {}
    pos_in_update = {}

    position = 0
    updatepos = 0


    in_use = []

    #Lists of 0s are constructed, which can then be updates to 1s if a given criteria is fulfilled
    #For example, if a gene is enriched in female flies, the 0 corresponding to the "female" category should be changed to 1
    #To this end, during default list construction, the position corresponding to each flytype is stored in the "Pos_in" dictionaries
    for flytype in weight_dict:
        
        #if the assigned weight for a given flytype is 0, it's one we're actually interested in. Thus, these contribute to the "true default"
        if weight_dict[flytype] > 0:
            in_use.append(flytype)
            true_default.append(0)
            pos_in_update[flytype] = updatepos
            updatepos += 1

        #The expression in other flytypes must also be taken into account, to allow "ALL" flytype to be analysed even if any other flytype is not desired
        pos_in_weights[flytype] = position
        position += 1
        default_gene_score.append(0)

    #The initial tally is created for each gene which has an A.gam ortholog
    #I recognise that defaultdict would be faster but would not allow for a flexible length of default (to my knowledge)
    for gene in has_an_ortho:
        moztubule_tally[gene] = [] + default_gene_score

    #The expression data for each gene is stored in the mtdata dictionary
    mtdata = extract_mt_expression(expressionfile)

    #as for MozAtlas, genes are then processed individually.
    for mozgene in moztubprobeids:
        fbgene = orthotodmel[mozgene]

        #Each probe is also handled seperately
        for probe in moztubprobeids[mozgene]:
            probedata = mtdata[probe]

            enrichment = {}
            
            #The whole fly reading is just stored as a variable
            wholefly = float(probedata[5])

            #For other fly type readings, the enrichment (avg expression / wholefly) is stored in the "enrichment" dictionary {Flytype:enrichment}
            enrichment['ADULT'] = float(probedata[2]) / wholefly
            enrichment['MALE'] = float(probedata[4])/wholefly
            enrichment['FEMALE'] = float(probedata[3])/wholefly
            enrichment['LARVAL'] = float(probedata[1])/wholefly

            #Each flytype then has its enrichment checked against the provided threshold, and successful fulfilment of criteria stored in moztubule_tally
            for flytype in weight_dict:
                
                #ALL must be handled a little differently - If EITHER the current probe enrichment > the threshold 
                # OR enrichment of any previous probe for any ortholog of the same D.mel gene is > the threshold, 
                #Then the moztubule_tally is updated to 1
                if flytype == 'ALL':
                    if (enrichment['ADULT'] > threshold or moztubule_tally[fbgene][pos_in_weights['ADULT']] == 1) and (enrichment['LARVAL'] > threshold or moztubule_tally[fbgene][pos_in_weights['LARVAL']] == 1):
                        moztubule_tally[fbgene][pos_in_weights[flytype]] = 1

                else:
                    if enrichment[flytype] > threshold:
                        moztubule_tally[fbgene][pos_in_weights[flytype]] = 1
            
        #Orthologs which are associated with multiple D.mel genes are accounted for
        if mozgene in problemorthos:
            for othergene in problemorthos[mozgene]:

                if othergene == fbgene:
                    pass

                else:
                    
                    #The tally for any other gene associated with a given ortholog is updated to the HIGHEST of the current tally or its existing tally.
                    #This applies for all flytypes
                    for flytype in weight_dict:
                        moztubule_tally[othergene][pos_in_weights[flytype]] = max([moztubule_tally[fbgene][pos_in_weights[flytype]], moztubule_tally[othergene][pos_in_weights[flytype]]])
                    
                    #Another check for the conditions for "ALL" being met is made (ADULT and LARVAL enrichment for any ortho of a given gene > threshold)
                    if 'ALL' in in_use and moztubule_tally[othergene][pos_in_weights['ADULT']] == 1 and moztubule_tally[othergene][pos_in_weights['LARVAL']] == 1:
                        moztubule_tally[othergene][pos_in_weights['ALL']] = 1
    
    #The desired readings are then transferred from the "moztubule_tally" to the "true_moztubule_tally" dictionary.
    #Handing results in this way allows us to check whether "ALL" flies are enriched > threshold
    for gene in has_an_ortho:
        updated_score = [] + true_default

        for flytype in in_use:
            updated_score[pos_in_update[flytype]] = moztubule_tally[gene][pos_in_weights[flytype]]

        true_moztubule_tally[gene] = updated_score
            
    return true_moztubule_tally, true_default, in_use

#This function prints MozAtlas results to the Mozatlas Breakdown .csv, for clarity to the user where a given tally comes from
def printoutscores(tally_compilation, dmeltoortho, threshold, threshold_spec, outfolder, allgenes):

    with open(f'{outfolder}/Agam_MozAtlas_Breakdown.csv', 'w+') as conservefile:
        OUTconserve = csv.writer(conservefile)
        
        #This block builds the header row
        header = ['FB_ID', 'Anopheles Gambiae ortholog(s)']

        for flytype in tally_compilation:
            header_extension = [f'Any Ortholog combination enriched >= {threshold} in target tissues vs whole mosquito - {flytype} flies', f'Any Ortholog combination at least {threshold_spec} times higher in target tissues compared to non target tissues? - {flytype} flies', f'ANY orthologs expressed in all target tissues? - {flytype} flies ']
            header.extend(header_extension)
        
        OUTconserve.writerow(header)

        #Gene data is printed: [D.mel gene, A.gam orthos (if any), then enrichment > threshold, abundance > non-target tissues, expression  in tissues of interest PER FLY TYPE]
        for gene in allgenes:

            try:
                if len(dmeltoortho[gene]) > 1:
                    agorthologs = '; '.join(dmeltoortho[gene])
                else:
                    agorthologs = dmeltoortho[gene][0]
            except:
                agorthologs = 'No orthologs found'

            outrow = [gene, agorthologs, ]

            for flytype in tally_compilation:
                outrow.extend(tally_compilation[flytype][gene])

            OUTconserve.writerow(outrow)

#This function prints MozTubule results to the Mozatlas Breakdown .csv, for clarity to the user where a given tally comes from      
def print_mt_outscores(moztubule_tally, dmeltoortho, threshold, outfolder, allgenes, mt_default, mt_in_use):

    with open(f'{outfolder}/Agam_MozTubules_Breakdown.csv', 'w+') as mtfile:
        OUTmt = csv.writer(mtfile)

        #This block builds the header row
        header = ['FB_ID', 'Anopheles Gambiae ortholog(s)',]

        for flytype in mt_in_use:
            header.append(f'Any combination of orthologs enriched >= {threshold} in {flytype} tubules?')

        OUTmt.writerow(header)

        #Gene data is printed: [D.mel gene, A.gam orthos (if any), then enrichment > threshold in tubules PER FLY TYPE]
        for gene in allgenes:

            try:
                if len(dmeltoortho[gene]) > 1:
                    agorthologs = '; '.join(dmeltoortho[gene])
                else:
                    agorthologs = dmeltoortho[gene][0]
            except:
                agorthologs = 'No orthologs found'

            outrow = [gene, agorthologs, ]

            try:
                outrow.extend(moztubule_tally[gene])
            except KeyError:
                outrow.extend(mt_default)

            OUTmt.writerow(outrow)

#This function uses the user-supplied weighting for the various readings (Enrichment, Abundance, Expression) and for MozTubules (compared to a tissue from MozAtlas)
#The final tally will be a score OUT OF ONE based on the entirity of Anopheles gambiae orthology
def agam_to_digittally(allgenes, moztubules_weight, tally_compilation, moztubules_tally, tissues_of_interest, measure_weights):

    #This will be a list of lists for csv printing
    gene_rows = []

    modified_weight = float(1)
    final_weight = float(1)

    #If both modules are desired, we need to prepare two weights (preset to 1 above)
    if moztubules_weight > 0 and moztubules_tally != {} and tally_compilation != {}:

        #Modified weight gets the weighting factor for the moztubule score based on the total number of tissues of interest
        #Ie, if only Tubules are input, by default, MozAtlas and MozTubules will be equally weighted
        modified_weight = float(moztubules_weight/len(tissues_of_interest))

        #The final weight is used to adjust the final score to the 0-1 scale
        final_weight =  float(1 + modified_weight)

    #The divisor_by_type generates a per-fly-type possible maximum score for MozAtlas
    divisor_by_type = (measure_weights['Enrichment'] + measure_weights['Specificity'] +  measure_weights['Expression'])

    for gene in allgenes:
        final_mozatlas_score = 0
        moztubules_score = 0
        
        #If MozAtlas has been investigated, we will adjust its scores based on user-defined component weights
        if tally_compilation != {}:
            total_divisor = 0
            mozatlas_score = 0

            #For each fly type, a score is compiled based on the weights of each ranking component.
            #Effectively, this score = the weight for each component if the requirements are met.
            #The total divisor = the divisor_by_type times the number of flytypes
            for flytype in tally_compilation:
                mozatlas_score += float(measure_weights['Enrichment'] * tally_compilation[flytype][gene][0])
                mozatlas_score += float(measure_weights['Specificity'] * tally_compilation[flytype][gene][1])
                mozatlas_score += float(measure_weights['Expression'] * tally_compilation[flytype][gene][2])

                total_divisor += divisor_by_type
            
            #The final mozatlas score = the total score across desired flytypes over the total divisor - ie, reduced to a 0-1 scale
            final_mozatlas_score = mozatlas_score/total_divisor

        #If MozAtlas has been investigated, we will adjust its scores based on user-defined MozTubule weight
        if moztubules_tally != {}:

            try:
                #If a gene is found in the moztubules_tally, a score is generated by dividing the number of desired flytypes which meet the enrichment threshold
                #by the total number of desired flytypes - ie, there is a max score of 1
                moztubules_score = float(sum(moztubules_tally[gene])/ len(moztubules_tally[gene]))

            except:
                #If the gene is not found in the tally, it has a score of 0 for all desired flytypes and as such a final score of 0
                moztubules_score = 0

            #The moztubules_score for each gene is then adjusted by the weight modified for MozTubules described above
            moztubules_score = moztubules_score * (modified_weight)

        #The FINAL score for each gene equals the score from MozAtlas (if any) plus the score from MozTubules) if any.
        #This is then divided by the final_weight parameter to adjust the scale to 0-1
        final_gene_score = float((final_mozatlas_score + moztubules_score)/final_weight)

        #A row to be printed for each gene is then stored in gene_rows
        gene_rows.append([gene, final_gene_score])
    
    return gene_rows

#This function creates the "Agam_orthology_tally csv, needed for the final DIGITtally output"
def append_to_tally(gene_rows, outfolder):
    with open(f'{outfolder}/Agam_orthology_tally.csv', 'w+') as tallyfile:

        OUTtally = csv.writer(tallyfile)

        header = ['Gene symbol', 'Anopheles gambiae score']
        OUTtally.writerow(header)

        for line in gene_rows:
            OUTtally.writerow(line)

def Score_Agam_Orthology(
        database, host, pw, user, orthologs, outputorthology, outputtally, moztubules_probes, moztubules_expression, tissues, atlasthresh_enr, 
        atlasthresh_spec, present_thresh, ma_targetweights, ma_measureweights, tubulesthresh, tub_to_atlas, mt_targetweights, modules
        ):
    
    #We read the user input weights into a dictionary (Key = Fly type, value = weight/total weights)
    mozatlas_weight_dict = gather_weights(ma_targetweights)
    mt_weight_dict = gather_weights(mt_targetweights)
    
    measure_weights = gather_weights(ma_measureweights)

    sexes_for_mozatlas = define_sexes(mozatlas_weight_dict)

    target_tissues = tissues

    orthologs, orthotodmel, dmeltoortho, problematicorthos, has_an_ortho, all_genes = pop_orthos(orthologs)

    tally_compilation = {}
    moztubule_tally = {}

    #The functions related to MozAtlas analysis are only run when MozAtlas is a desired function
    if 'MozAtlas' in modules:
        cnx = create_server_connection(host, user, pw, database)

        probeids = get_probe_ids(cnx, orthologs)

        categories = find_mozatlas_expression(cnx, has_an_ortho, sexes_for_mozatlas, target_tissues, probeids, orthotodmel, atlasthresh_enr, atlasthresh_spec, present_thresh, problematicorthos)

        tally_compilation = score_by_cat(categories, mozatlas_weight_dict)

        printoutscores(tally_compilation, dmeltoortho, atlasthresh_enr, atlasthresh_spec, outputorthology, all_genes)

    #The functions related to MozTubule analysis are only run when Malpighian Tubules are a tissue of interest and MozTubules is a desired function
    if 'MozTubules' in modules:
        print('MozTubules starting')
        moztubprobeids = get_moztubule_probeids(orthologs, moztubules_probes)

        moztubule_tally, mt_default, mt_in_use = find_moztubules_expression(moztubprobeids, orthotodmel, tubulesthresh, moztubules_expression, problematicorthos, mt_weight_dict, has_an_ortho)
        
        print_mt_outscores(moztubule_tally, dmeltoortho, tubulesthresh, outputorthology, all_genes, mt_default, mt_in_use)

        print('Moztubules Complete')
        
    #The DIGITtally output is ALWAYS produced
    rows_to_write = agam_to_digittally(all_genes, tub_to_atlas, tally_compilation, moztubule_tally, target_tissues, measure_weights)
    append_to_tally(rows_to_write, outputtally)
