#Code written by Dr. Andrew Gillen (Dow-Davies lab, University of Glasgow)
#Originally written October 2022, last updated 26/06/2023

import csv
import os
import glob

from collections import defaultdict

#Part of the pipeline for DIGITtally analysis - see www.digittally.org
#Combines data from FlyAtlas1 and FlyAtlas2 into a single file which will serve as the basis for DIGITtally scoring
#Importantly, this also creates a list of genes which pass one of the user-defined thresholds in either FA1 or FA2 - These are the genes of interest for this run

#A simple function to break up the weight string (-tw arg) supplied into an accessible dictionary
def gather_weights(targetweights):
    weight_dict_untouched = {}
    weight_dict = {}
    weight_dict['FlyAtlas1'] = {}
    weight_dict['FlyAtlas2'] = {}
    obligates = []
    weights = targetweights.split(';')

    sumweightsfa1 = 0
    sumweightsfa2 = 0
    adultfa1done = 0

    #parses the supplied weight into a dictionary
    for item in weights:
        subitem = item.split(':')
        indivweight = float(subitem[1])
        weight_dict_untouched[subitem[0]] = indivweight
        sumweightsfa2 += indivweight
        
        if subitem[0] in ['MALE', 'FEMALE', 'ADULTS']:
            sumweightsfa1 += indivweight
        else:
            sumweightsfa1 += indivweight
        
        #This facilitates designating a fly type as "Obligatory" - ie, a gene MUST be expressed in this type to be listed as a gene of interest
        if len(subitem) == 3:
            if subitem[2] == 'ob':
                obligates.append(subitem[0])
    
    #Divides each weight by the sum of all weights provided, allowing fine tuning of weighting of different fly types.
    for flytype in weight_dict_untouched:
        weight_dict['FlyAtlas2'][flytype] = (weight_dict_untouched[flytype]/sumweightsfa2)

        if flytype in ['MALE', 'FEMALE', 'ADULTS']:
            if 'ADULT' in  weight_dict['FlyAtlas1']:
                weight_dict['FlyAtlas1']['ADULT'] = max([weight_dict['FlyAtlas1']['ADULT'],(weight_dict_untouched[flytype]/sumweightsfa1)])
            else:
                weight_dict['FlyAtlas1']['ADULT'] = (weight_dict_untouched[flytype]/sumweightsfa1)
        else:
            weight_dict['FlyAtlas1'][flytype] = (weight_dict_untouched[flytype]/sumweightsfa1)

    return weight_dict, obligates

#using the FlyBase FBgn ⇔ FBtr ⇔ FBpp IDs (expanded) file, the information for each gene (Symbol, transcripts) are populated to the appropriate FbID
def populate_gene_info(file):
    genebytranscript = {}
    symbolbygene = {}

    with open(file) as transcriptfile:
        INtranscripts = csv.reader(transcriptfile, delimiter='\t')

        for line in INtranscripts:
            
            #Comment lines are skipped over, the rest are parsed to provide symbols and transcripts by gene id
            if line[0][0] != '#':
                geneid = line[2]
                transcriptid = line[7]
                symbol = line[3]

                genebytranscript[transcriptid] = geneid

                if geneid not in symbolbygene:
                    symbolbygene[geneid] = symbol

    return genebytranscript, symbolbygene

#Taking the inputs from FlyAtlas Analyses, a "hitlist" as created containing all GENE IDs associated with any of the user-defined lists of interest
def get_hits(hitlist, file):
    specificcat = []

    with open(file) as processing:
        for line in processing:

            possiblehit = line.strip('\n')

            #Sometimes, a FlyAtlas1 reading cannot be definitively assigned to a single gene. In these cases, all possible genes are added to the hit list
            if '///' in possiblehit:
                partialpos = possiblehit.split(' /// ')
                
                for partial in partialpos:
                    if partial not in hitlist:
                        hitlist.append(partial)

                    specificcat.append(partial)

            else:        
                if possiblehit not in hitlist:
                    hitlist.append(possiblehit)
                
                specificcat.append(possiblehit)
    
    return hitlist, specificcat

#Functions as above, but is specifically designed to handle FlyBase transcript IDs, converting them into Gene IDs if possible and appending them to a "not found" list if not.
def transcript_hits(hitlist, file, genebytranscript, specifictranscripts, needsasynonym):
    specificcat = []

    with open(file) as processing:
        for line in processing:
            
            hittranscript = line.strip('\n')

            #We use the dictionary generated earlier to find gene by transcript ID
            try:
                hitgene = genebytranscript[hittranscript]

                if hitgene not in hitlist:
                    hitlist.append(hitgene)
                
                specificcat.append(hitgene)

                if hitgene not in specifictranscripts:
                    specifictranscripts[hitgene] = set()

                specifictranscripts[hitgene].add(hittranscript)
            
            #If the transcript ID is not found in the dictionary, it is out of date and a more recent identifier must be found to allow further study
            except:
                needsasynonym['T'].add(hittranscript)
    
    return hitlist, specificcat, specifictranscripts, needsasynonym

#A function to check for alternate identifiers for genes which cannot be found
def get_synonyms(needsasynonym, hitlist, synonymfile, outfolder):

    hitlist_modified = [] + hitlist

    genestosearch = needsasynonym['G']
    transcriptstosearch = needsasynonym['T']
    updateids = {}
    old_to_new = {}
    found = []

    #The FlyBase_SecondarySynonyms.tsv is used to find alternative IDs for each gene, which can be recognised by contemporary FlyBase.
    #It is unclear for how long a deprecated ID is stored by FlyBase, and so even with this measure some genes may not be found
    with open(synonymfile) as synfile, open(f'{outfolder}/Synonyms/DetectedGOIs_FlyBaseIDUpdates.csv', 'w+') as detectfile:
        INsyn = csv.reader(synfile, delimiter = '\t')
        OUTdetected = csv.writer(detectfile)

        header = ['FbID in FlyAtlas1/2', 'UPDATED ID from FlyBase']
        OUTdetected.writerow(header)

        for line in INsyn:

            if len(line) > 0:

                if line[0][0] != '#': 

                    newid = line[2]

                    for igene in genestosearch:

                        if igene in line[3]:

                            if newid not in hitlist_modified:

                                hitlist_modified.append(newid)
                            
                            old_to_new[igene] = [newid]
                            updateids[newid] = igene
                            found.append(igene)

    #The IDs of unfindable genes are printed
    with open(f'{outfolder}/Unfindable_FbIDs.txt', 'w+') as missgenefile:
        nunf = 0

        for gene in genestosearch:

            if gene not in found:
                #As FlyBase continues to update, if FlyAtlases are not kept up to date, this will be hit more and more frequently
                missgenefile.write(f'{gene}\n')
                nunf += 1

            while gene in hitlist_modified:
                hitlist_modified.remove(gene)
    
    if nunf > 0:
        print(f'\nEven with the secondary symbols, {nunf} cannot be satisfactorally married up to a current FlyBase ID. As such, these will need to be excluded from DIGITtally, though may prompt further research')
        print(f'Unfindable FbIDs will be printed to {outfolder}/Unfindable_FbIDs.txt to allow user-directed study')          
    else:
        os.remove(f'{outfolder}/Unfindable_FbIDs.txt')

    #Transcripts cannot currently be updated to more recent versions as they are not stored by FlyBase. Ways to work around this are being investigated
    if len(transcriptstosearch) > 0:
        print('\nUnfortunately, there is currently no way to marry up out-of-date transcript IDs with current identifiers')
        print(f'For the moment, Transcripts which cannot be married up to a current FlyBase ID will be printed to {outfolder}/Out-of-Date_TranscriptHits.txt, allowing user-directed study')

        #Unfound IDs are stored in a seperate file
        with open(f'{outfolder}/Out-of-Date_TranscriptHits.txt', 'w+') as missfile:
            for item in transcriptstosearch:
                missfile.write(f'{item}\n')

    return hitlist_modified, updateids, old_to_new

#This function checks which genes from a FlyAtlas analysis pass user-defined flytype obligations
def check_against_needs(flytype_sub_dictionary, weight_dict_for_source, outputbyhit, tmpobs, obstarget, updateids):

    currentscore = {}

    obligationsmet = {}

    for flytype in flytype_sub_dictionary:

        if weight_dict_for_source[flytype] != 0:

            for hit in outputbyhit:
                datafound = 0

                if hit not in currentscore:
                    currentscore[hit] = 0
                    obligationsmet[hit] = 0

                weightmod = weight_dict_for_source[flytype]

                if hit in flytype_sub_dictionary[flytype]:
                    currentscore[hit] += (1 * weightmod)
                    datafound = 1

                    #This flags up when a gene is expressed in an obligatory fly type
                    if flytype in tmpobs:
                        obligationsmet[hit] += 1

                #This catches edge cases where the old ID is used in FlyAtlas1 but the newer ID in FlyAtlas2
                if hit in updateids and datafound == 0:
                    if updateids[hit] in flytype_sub_dictionary[flytype]:
                        currentscore[hit] += (1 * weight_dict_for_source[flytype])
                        datafound = 1
                    
                        if flytype in tmpobs:
                            obligationsmet[hit] += 1

    #This block checks that every putative gene of interest passes the measurement metric in all obligate fly types 
    #If not, the score for that metric is reduced to 0
    for hit in currentscore:
        if obligationsmet[hit] >= obstarget:
            outputbyhit[hit].append(currentscore[hit])
        else:
            outputbyhit[hit].append(0)
    
    return outputbyhit

#This function created the "Tally Starter", which will be used as the base for the eventual full DIGITtally output
def build_tally_starter(hitlist, tallycats, symbols, specifictranscripts, updateids, max_len, weight_dict, obligates, completed_thresholds, header, additional_threshold_data, outfolder, user_supplied_list):
    
    outputbyhit = {}

    mode_to_label = { 'ENRICHED' : 'ENRICHMENT', 'ABUNDANT' : 'ABUNDANCE'}

    #Ensures each hit (flybase ID) is joined to a symbol and, where appropriate, the specific transcripts enriched in tissues of interest
    for hit in hitlist:
        outputbyhit[hit] = [hit, symbols[hit]]

        try:
            outputbyhit[hit].append('; '.join(specifictranscripts[hit]))
        except:
            outputbyhit[hit].append('n/a')
        
        outputbyhit[hit].append(' ')

    #Each data source which is to be included is used to generate seperate scores
    for datasource in tallycats:
        tmpobs = []
        tmpobs += obligates[datasource]

        #The list of obligates must also be adapted for the ADULT-only nature of FlyAtlas1
        if datasource == 'FlyAtlas1':
            obstarget = sum([1 for x in obligates[datasource] if x not in ['MALE', 'FEMALE', 'ADULTS']])

            addition = 0

            for flytype in ['MALE', 'FEMALE', 'ADULTS']:
                if flytype in obligates[datasource]:
                    addition = 1
                    tmpobs.append('ADULT')
            
            obstarget += addition

        else:
            obstarget = len(obligates[datasource])
        
        #scores are generated for each desired metric (enrichment and/or abundance)
        for mode in tallycats[datasource]:
            header.append(f'{datasource}_{mode}_{completed_thresholds[datasource][mode_to_label[mode]][0]}_SCORE')

            #These scores are generated using a tally weighted based on the user defined weight
            #Score = sum of all fly types where individual fly types scored as (1 if detected in a tissue) * weight for fly type
            #We need to explicitly set aside 'Adult' FA1 types so the weights can be processed correctly. 
            outputbyhit = check_against_needs(tallycats[datasource][mode], weight_dict[datasource], outputbyhit, tmpobs, obstarget, updateids)

        if datasource in additional_threshold_data:

            for mode in additional_threshold_data[datasource]:

                for extra_threshold in additional_threshold_data[datasource][mode]:
                    header.append(f'{datasource}_{mode}_{extra_threshold}_SCORE')

                    outputbyhit = check_against_needs(additional_threshold_data[datasource][mode][extra_threshold], weight_dict[datasource], outputbyhit, tmpobs, obstarget, updateids)

    #This catches cases (generally broad investigations eg of ONE tissue) where > the max number of desired genes of interest are found.
    #An upper limit on genes of interest is generally desirable as large gene lists will greatly increase execution times of later programs
    #If you're willing to wait, this can be circumvented by supplying an arbitrarily large (>20000) -max argument
    if len(outputbyhit) > max_len:
        scorebygene = []
        hitlistupdate = []

        for gene in outputbyhit:
            scorebygene.append((gene, sum(outputbyhit[gene][4:])))
        
        scorebygene.sort(key=lambda i:i[1],reverse=True)
        
        usethese = [x[0] for x in scorebygene[:max_len]]

        for gene in outputbyhit:
            if gene in usethese:
                hitlistupdate.append(gene)
        
        print(f'Note that, due to large numbers of hit genes and program settings, only the {max_len} best "scoring" genes of these hits are included in your DIGITtally \n\n\n')

    else:
        hitlistupdate = hitlist

    #We remove all genes which have a sum of all scores of 0 - these are genes which WOULD be genes of interest in the defined tissues but DO NOT meet the fly type obligations
    update2 = []
    for item in hitlistupdate:
        if sum(outputbyhit[item][4:]) > 0 or user_supplied_list == 1:
            update2.append(item)
    
    hitlistupdate = [] + update2
    hitlistupdate = set(hitlistupdate)

    #Tallies are printed to the DIGITtally starter file for the top x genes of interest.
    with open(f'{outfolder}/DIGITtallyStarter.csv', 'w+') as tallyfile:
        OUTtally = csv.writer(tallyfile)
        
        OUTtally.writerow(header)

        for hit in hitlistupdate:
            OUTtally.writerow(outputbyhit[hit])

    return hitlistupdate, outputbyhit, header

#simple .txt file lists of gene of interest printed for use in later programming
def make_lists(hitlist, symbols, outfolder):
    with open(f'{outfolder}/DetectedGOIs_FbID.txt', 'w+') as fbidfile, open(f'{outfolder}/DetectedGOIs_Symbol.txt', 'w+') as symbolfile:
        
        for hitgene in hitlist:
            hitsymbol = symbols[hitgene]

            fbidfile.write(f'{hitgene}\n')
            symbolfile.write(f'{hitsymbol}\n')

#This is used for extra threshold processing - filenames are used to inform the agglomerator on thresholds
def update_thresholds_from_filename(filename, dictionary):
    parsed = filename.split('_')

    abn_thresh = parsed[-1]
    abn_thresh = abn_thresh.replace('-', '.')
    enr_thresh = parsed[-2]
    enr_thresh = enr_thresh.replace('-', '.')

    dictionary['ABUNDANCE'].append(abn_thresh)
    dictionary['ENRICHMENT'].append(enr_thresh)

    return dictionary

#This builds a dictionary which stores lists of abundance/enrichment thresholds
def initialise_threshold_checking(fa1, fa2):
    completed_thresh_dict = {}

    for atlas in ['FlyAtlas1', 'FlyAtlas2']:
        completed_thresh_dict[atlas] = defaultdict(list)

    try:
        completed_thresh_dict['FlyAtlas1'] = update_thresholds_from_filename(fa1, completed_thresh_dict['FlyAtlas1'])
    except:
        pass

    try:
        completed_thresh_dict['FlyAtlas2'] = update_thresholds_from_filename(fa2, completed_thresh_dict['FlyAtlas2'])
    except:
        pass

    return completed_thresh_dict

#This function scores genes which pass non-base thresholds. As these thresholds are non-base, they ARE NOT used to inform hitlist construction (as they would not be informative.)
def process_extra_thresholds(folder_format, atlas, completed_thresholds, desigmodes, targets, weight_dict, genebytranscript):

    #{data type: {threshold: {fb_id:score}}}
    extra_data = defaultdict(lambda: defaultdict(dict))
    mode_to_label = { 'ENRICHED' : 'ENRICHMENT', 'ABUNDANT' : 'ABUNDANCE'}

    for folder in glob.glob(f'{folder_format}*'):

        new_thresholds = defaultdict(list)
        new_thresholds = update_thresholds_from_filename(folder, new_thresholds)
        
        process_threshold = defaultdict(lambda:0)

        for data_type in desigmodes:

            for thresh_value in new_thresholds[mode_to_label[data_type]]:

                if thresh_value not in completed_thresholds[atlas][mode_to_label[data_type]]:
                    process_threshold[data_type] = 1

        temp_data_fragment = defaultdict(dict)

        if atlas == 'FlyAtlas1':
            unused_hitlist, temp_data_fragment = handle_fa1_data(targets, weight_dict, desigmodes, folder)

        elif atlas == 'FlyAtlas2':
            unused_hitlist, temp_data_fragment, unused_needsasynonym, unused_specifictranscripts = handle_fa2_data(targets, weight_dict, desigmodes, genebytranscript, folder)


        for data_type in temp_data_fragment:

            if process_threshold[data_type] == 1:
                extra_data[data_type][new_thresholds[mode_to_label[data_type]][0]] = temp_data_fragment[data_type]
    
    return extra_data
      
#Gathers genes of interest and their scores from FlyAtlas1  
def handle_fa1_data(targets, weight_dict, desigmodes, folder):

    #This exists to catch a VERY SPECIFIC issue - as FlyAtlas1 only contains generic "Adult" data, feeding in weights for male or female defaults to the adult reading
    #To prevent multiple reading of the same data artificially increasing the weight given to FlyAtlas1 Adults, these cases are only to be counted ONCE
    exceptionadultadded = 0
    
    hitlist = []

    datadict = defaultdict(dict)

    #For every fly type, gene of interest lists from FlyAtlas1_Analyser will be run
    for flytype in targets:
        #Only fly types we're interested in, ie types with a weight > 0, will be processed.
        try:
            if weight_dict['FlyAtlas1'][flytype] != 0:
                
                flytypetosearch = flytype

                #Genes of interest which are Enriched or Abundant (as desired) will be processed seperately, providing two seperate scores for each data source.
                for mode in desigmodes:

                    fa1filetocheck = f'{folder}/{flytypetosearch}_{mode}.txt'

                    hitlist, specificcat = get_hits(hitlist, fa1filetocheck)
                    
                    datadict[mode][flytype] = specificcat

        #As noted above - if MALE OR FEMALE are being investigated, we default to the ADULTS genes of interest while storing the fact this has happened
        except:
            if (flytype in ['MALE', 'FEMALE', 'ADULTS'] and weight_dict['FlyAtlas1']['ADULT'] != 0):
                flytypetosearch = 'ADULTS'
                exceptionadultadded += 1

                for mode in desigmodes:

                    #If the ADULTS data has been processed once, it will not be saved again
                    if exceptionadultadded < 2:

                        fa1filetocheck = f'{folder}/{flytypetosearch}_{mode}.txt'

                        hitlist, specificcat = get_hits(hitlist, fa1filetocheck)

                        datadict[mode]['ADULT'] = specificcat

    return hitlist, datadict

#Gathers genes of interest and their scores from FlyAtlas2
def handle_fa2_data(targets, weight_dict, desigmodes, genebytranscript, folder, needsasynonym = {'T':set()}):
    
    hitlist = []

    specifictranscripts = {}

    datadict = defaultdict(dict)

    for flytype in targets:

        if weight_dict['FlyAtlas2'][flytype] != 0:

            for mode in desigmodes:

                fa2genefile = f'{folder}/FPKMGene_{mode}in{flytype}.txt'

                hitlist, specificcatgenes = get_hits(hitlist, fa2genefile)
                
                fa2transcriptfile = f'{folder}/FPKMTranscript_{mode}in{flytype}.txt'
                hitlist, specificcattranscripts, specifictranscripts, needsasynonym = transcript_hits(hitlist, fa2transcriptfile, genebytranscript, specifictranscripts, needsasynonym)

                #We compile genes detected at the "whole gene" level and at the "transcript" level into ONE list - 
                #these are not added multiple times if eg the whole gene is tissue-enriched and ALSO a specific transcript is tissue-enriched
                
                specificcatcompiled = []
                specificcatcompiled = specificcatcompiled + specificcatgenes

                for gene in specificcattranscripts:

                    if gene not in specificcatcompiled:

                        specificcatcompiled.append(gene)
                
                
                datadict[mode][flytype] = specificcatcompiled

    return hitlist, datadict, needsasynonym, specifictranscripts

#This builds weight dicts containing all 0s for cases where the user has uploaded their own gene list and excluded both FlyAtlases
def generate_default_weights():
    default_weight_dict = {
        'FlyAtlas1':defaultdict(lambda:0),
        'FlyAtlas2':defaultdict(lambda:0)
    }

    return default_weight_dict

def make_tally_basics(preexisting_list, fa1_folder, fa2_folder, fa1extra, fa2extra, mode, outputlists, outputtally, infofile, secondaryannotations, weightsfa1, weightsfa2, maxlen):
    print('STARTING')
    weight_dict = generate_default_weights()
    print(weight_dict)
 
    #Some definitions needed
    final_hitlist = []
    user_supplied_list = 0
    targets = ['MALE', 'FEMALE', 'LARVAL', 'ADULTS', 'ALL']

    tallycats = defaultdict(lambda: defaultdict(dict))
    obligates = {}

    specifictranscripts = {}

    completed_thresholds = initialise_threshold_checking(fa1_folder, fa2_folder)

    needsasynonym = {}
    for letter in ['T', 'G']:
        needsasynonym[letter] = set()

    #We gather info on the transcripts and Gene Symbol ssociated with each gene identifier
    genebytranscript, symbolbygene = populate_gene_info(infofile)

    #Sets up a list of hit types to investigate.
    if mode == 'both':
        desigmode = ['ENRICHED', 'ABUNDANT']
        skipparameter = False
        
    elif mode == 'neither':
        skipparameter = True

    else:
        desigmode = [f'{mode.upper()}']
        skipparameter = False


    #...If a FlyAtlas1 folder is provided, allowing the program to be skipped if not desired
    if fa1_folder != '' and skipparameter == False:

        #We read the user input weights into a dictionary (Key = Fly type, value = weight/total weights)
        weight_dict, obligates['FlyAtlas1'] = gather_weights(weightsfa1)

        fa1hits, tallycats['FlyAtlas1'] = handle_fa1_data(targets, weight_dict, desigmode, fa1_folder)

        final_hitlist.extend(fa1hits)

    else:
        obligates['FlyAtlas1'] = []

    #A very similar block for FlyAtlas2_Analyser files. Works the same as above aside from two differences:
    #No need to explicitly prevent data re-processing
    #Transcripts of interest must be processed.
    if fa2_folder != '' and skipparameter == False:
        #We read the user input weights into a dictionary (Key = Fly type, value = weight/total weights)
        weight_dict, obligates['FlyAtlas2'] = gather_weights(weightsfa2)

        fa2hits, tallycats['FlyAtlas2'], needsasynonym, specifictranscripts = handle_fa2_data(targets, weight_dict, desigmode, genebytranscript, fa2_folder, needsasynonym)

        final_hitlist.extend(fa2hits)
    else:
        obligates['FlyAtlas2'] = []

    if preexisting_list != []:
        user_supplied_list = 1
        final_hitlist = preexisting_list

    #We check to see if we can associate each gene of interest FbID with a Gene Symbol. 
    #If not, we need to search for an alternate ID number for a symbol
    for hit in final_hitlist:
        try:
            temp = symbolbygene[hit]

        except:
            needsasynonym['G'].add(hit)

    #Here we find backup IDs for genes, and use these to correct our hitlist to the most up to date FlyBase identifiers. This is necessary for querying FlyBase data later.
    #THEORETICALLY we would also do this for transcripts which cannot be identified using their FlyAtlas2 identifiers - however, FlyBase does not seem to keep track of previous transcript IDs
    final_hitlist, updateids, old_to_new = get_synonyms(needsasynonym, final_hitlist, secondaryannotations, outputlists)
    
    header = ['FbID', 'Gene Symbol', 'Specific transcript(s) of interest','  ']
    
    additional_threshold_data = {}

    if fa1extra != '':
        additional_threshold_data['FlyAtlas1'] = process_extra_thresholds(fa1extra, 'FlyAtlas1', completed_thresholds, desigmode, targets, weight_dict, genebytranscript)
    
    if fa2extra != '':
        additional_threshold_data['FlyAtlas2'] = process_extra_thresholds(fa2extra, 'FlyAtlas2', completed_thresholds, desigmode, targets, weight_dict, genebytranscript)

    #Using the data we've gathered, we create output files
    try:
        final_hitlist, output_by_hit, header = build_tally_starter(final_hitlist, tallycats, symbolbygene, specifictranscripts, updateids, maxlen, weight_dict, obligates, completed_thresholds, header, additional_threshold_data, outputtally, user_supplied_list)
        make_lists(final_hitlist, symbolbygene, outputlists)
    except Exception as e:
        print(e)
