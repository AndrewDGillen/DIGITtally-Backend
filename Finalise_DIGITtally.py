import csv
import datetime

from collections import defaultdict

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

#Generates dictionaries contaning: weighting for BROAD DIGITtally categories; and for SPECIFIC data sources
def populate_weight_dicts(broad, specific):

    broad_weights = gather_weights(broad)

    specific_weights = {}
    pass_source = []

    for data_source in specific:
        specific_weights[data_source[0]] = gather_weights(data_source[1])

        if specific_weights[data_source[0]]['GLOBAL'] == 0:
            pass_source.append(data_source[0])

    return broad_weights, specific_weights, pass_source

#Populates information from the "DIGITtally starter" file generated from GOIAgglomerator.
#IE FA1/FA2 data, enrichment and abundance
def get_started(folder):

    data_by_id = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

    all_hit_ids = []

    id_by_symbol = {}
    symbol_by_id = {}
    transcripts_by_id = {}

    with open(f'{folder}/DIGITtallyStarter.csv') as starterfile:

        INstarter = csv.reader(starterfile)

        line_count = 0

        info_by_pos = defaultdict(lambda: ['', '', ''])

        for line in INstarter:
            position = 0

            if line_count == 0:

                for item in line:
                    
                    if position > 3:
                        split_info = item.split('_')

                        info_by_pos[position][1] = split_info[0]
                        info_by_pos[position][2] = split_info[2]

                        if split_info[1] == "ENRICHED":
                            info_by_pos[position][0] = 'Enrichment'
                        elif split_info[1] == "ABUNDANT":
                            info_by_pos[position][0] = 'Specificity'

                    position += 1

            #Given this file contains a full list of our genes of interest, we use it to populate a couple of things
                #- A list of IDs by gene symbol
                #- A list of all hit IDs
                #- A list of Gene Symbols By ID
                #- A list of transcripts of interest by gene ID
            else:
                fb_id = line[0]
                symbol = line[1]

                id_by_symbol[symbol] = fb_id
                
                all_hit_ids.append(fb_id)

                symbol_by_id[fb_id] = symbol
                transcripts_by_id[fb_id] = line[2]

                for item in line:

                    if position > 3:
                        data_type = '' + info_by_pos[position][0] + '_' + str(info_by_pos[position][2])
                        source = '' + str(info_by_pos[position][1]) 
                        data_by_id[data_type][source][fb_id] = item

                    position += 1                     

            line_count += 1

    return data_by_id, id_by_symbol, symbol_by_id, transcripts_by_id, all_hit_ids

#Populates all possible synonyms for a gene, can be necessary when converting formats for SCOPE etc
def get_synonyms(synonym_input):
    
    synonyms = defaultdict(list)

    try:
        with open(synonym_input) as synonym_file:
            INsynonyms = csv.reader(synonym_file)

            next(INsynonyms)

            for line in INsynonyms:

                synonyms[line[1]].append(line[0])

    except:
        print('No appropriate synonym file supplied')
    
    return synonyms

#Converts a gene symbol to a flybase ID
def symbol_to_id(symbol, id_by_symbol, synonyms):
    try:
        fb_id = [id_by_symbol[symbol]]

    #If the gene symbol cannot be found, we check the known synonym file
    except:

        try:
            synonym = synonyms[symbol]
            fb_id = [id_by_symbol[i_synonym] for i_synonym in synonym]
        except:
            print(f'{symbol} cannot be recognised by the supplied synonyms. Are you sure you supplied the correct synonym file?')
            exit()
    
    return fb_id

#Populates the data based on FlyCellAtlas findings - 5 different potential data streams. ALL of these are handled, then weighted according to user wants
def expand_fca(folder, data_by_id, id_by_symbol, synonyms, all_hit_ids):
    
    genes_in_fca = set()

    with open(f'{folder}/FlyCellAtlas_tally.csv') as fcafile:
        INfca = csv.reader(fcafile)

        next(INfca)

        fca_data_types = [('Specificity_NonTarget Expression', 2), ('Specificity_Proportion Expressing Gene Are Target', 3),('Ubiquity', 1),('Enrichment', 4)]

        for line in INfca:
            symbol = line[0]

            fb_id = symbol_to_id(symbol, id_by_symbol, synonyms)

            for data_type in fca_data_types:

                for indiv_id in fb_id:
                    
                    genes_in_fca.add(indiv_id)

                    data_by_id[data_type[0]]['FlyCellAtlas'][indiv_id] = line[data_type[1]]
            
            try:
                for indiv_id in fb_id:
                    data_by_id['Coexpression']['FlyCellAtlas'][indiv_id] = line[5]
            except:
                for indiv_id in fb_id:
                    data_by_id['Coexpression']['FlyCellAtlas'][indiv_id] = 0
    
    #deals with cases which cannot be searched in SCOPE, in which case genes score 0 for all measures
    for indiv_id in all_hit_ids:

        if indiv_id not in genes_in_fca:

            for data_type in fca_data_types:
                data_by_id[data_type[0]]['FlyCellAtlas'][indiv_id] = 0

            data_by_id['Coexpression']['FlyCellAtlas'][indiv_id] = 0
    
    return data_by_id

#Populates FlyBase information from the FlyBase tally
def expand_flybase(folder, data_by_id):

    with open(f'{folder}/FlyBase_tally.csv') as fbfile:

        INfb = csv.reader(fbfile)
        
        next(INfb)

        for line in INfb:

            fb_id = line[0]

            data_by_id['Known Activity_Association']['FlyBase'][fb_id] = line[1]
            data_by_id['Known Activity_Phenotype']['FlyBase'][fb_id] = line[2]
    
    return data_by_id

#Populates ALL orthology information, as these are BROADLY in a similar format
def process_orthology(folder, data_by_id, species, id_by_symbol, synonyms):

    species_parsed = species.split(' ')
    species_short = (species_parsed[0][0]) + (species_parsed[1][0:3])
    
    with open(f'{folder}/{species_short}_orthology_tally.csv') as orthologyfile:
        INortho = csv.reader(orthologyfile)

        next(INortho)

        #All species use a standardised format for their tally files, so we can easily grab the number we need.
        for line in INortho:
                fb_id = line[0] 
                data_by_id[f'Orthology_{species_short}']['ortho'][fb_id] = line[1]


    return data_by_id

#Now that results dicts are populated, we weight based on what modules the user ACTUALLY wants to contribute to their score
def weight_results(data_by_id, broad_weights, specific_weights):
    weighted_data = defaultdict(list)
    print(broad_weights)
    header = ['FlyBase ID', 'Gene Symbol', 'Specific Transcripts found in tissues of interest (from FlyAtlas2)', ' ']

    score_by_id = defaultdict(lambda: 0)
    list_of_scores = []
    
    to_sort = [x for x in data_by_id]

    sorted_list = sorted(to_sort)

    for category in sorted_list:

        if '_' in category:
            category_split = category.split('_')

            cat_heading = category_split[0]
            supplemental_info = category_split[1]
        
        else:
            cat_heading = category
            supplemental_info = ''

        broad_cat_weight = broad_weights[cat_heading]

        if broad_cat_weight != 0:
            for data_source in data_by_id[category]:

                global_source_weight = specific_weights[data_source]['GLOBAL']
                
                if cat_heading == 'Orthology':

                    source_cat_weight = specific_weights['ortho'][supplemental_info]
                    
                elif cat_heading == 'Known Activity':
                    source_cat_weight = specific_weights['FlyBase'][supplemental_info]
                    
                else:
                    source_cat_weight = specific_weights[data_source][cat_heading]

                if global_source_weight != 0 and source_cat_weight != 0:

                    if cat_heading != 'Orthology':
                        data_description = '' + cat_heading + '_' + data_source

                        if supplemental_info != '':
                            data_description = data_description + '_' + supplemental_info

                    else:
                        data_description = '' + supplemental_info + '_' + category
                        
                    header.append(data_description)

                    for indiv_id in data_by_id[category][data_source]:
                        
                        tally = float(data_by_id[category][data_source][indiv_id])
                        
                        weighted_tally = tally * broad_cat_weight * global_source_weight * source_cat_weight

                        weighted_data[indiv_id].append(weighted_tally)

                        score_by_id[indiv_id] += weighted_tally
    
    #Having build score by ID, we create a list of (ID, SCORE tuples)
    for fb_id in score_by_id:
        list_of_scores.append((fb_id, score_by_id[fb_id]))
    
    #Then sort the tuples to produce an ordered list
    ordered_scores = sorted(list_of_scores, key = lambda tup:tup[1], reverse = True)
    return weighted_data, header, ordered_scores

#For each gene, we gather up all the information we have and print it out in a single row, ordered by total tally score
def create_digittally_output(folder, weighted_data, tally_header, symbol_by_id, transcripts_by_id, order):

    with open(f'{folder}/DIGITtally_FINAL.csv', 'w+') as finalfile:
        OUTdigittally = csv.writer(finalfile)

        digitally_length = len(tally_header) - 4
        tally_header.append(' ')
        tally_header.append(f'DIGITtally (out of {digitally_length} MAXIMUM)')

        OUTdigittally.writerow(tally_header)

        for individual_gene in order:
            fb_id = individual_gene[0]
            digittally = individual_gene[1]

            outrow = [fb_id, symbol_by_id[fb_id], transcripts_by_id[fb_id], ' ']
            outrow.extend(weighted_data[fb_id])

            outrow.append(' ')
            outrow.append(digittally)

            OUTdigittally.writerow(outrow)

def WrapUp(tally_file_folder, synonym_file, broad_categories, userupload, flyatlas1, flyatlas2, flycellatlas, flybase, orthology, species):

    specific_data_sources = [("UserUpload", userupload), ('FlyAtlas1', flyatlas1), ('FlyAtlas2', flyatlas2), ('FlyCellAtlas', flycellatlas), ('FlyBase', flybase), ('ortho', orthology)]

    broad_weights, specific_weights, passsource = populate_weight_dicts(broad_categories, specific_data_sources)
    
    data_by_id, id_by_symbol, symbol_by_id, transcripts_by_id, all_hit_ids = get_started(tally_file_folder)

    if synonym_file != '':
        synonyms = get_synonyms(synonym_file)

    if 'FlyCellAtlas' not in passsource:
        data_by_id = expand_fca(tally_file_folder, data_by_id, id_by_symbol, synonyms, all_hit_ids)

    if 'FlyBase' not in passsource:
        data_by_id = expand_flybase(tally_file_folder, data_by_id)

    for indiv_species in species:

        data_by_id = process_orthology(tally_file_folder, data_by_id, indiv_species, id_by_symbol, synonyms)
    
    weighted_data, tally_header, order = weight_results(data_by_id, broad_weights, specific_weights)

    create_digittally_output(tally_file_folder, weighted_data, tally_header, symbol_by_id, transcripts_by_id, order)
