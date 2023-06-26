#Code written by Dr. Andrew Gillen (Dow-Davies lab, University of Glasgow)
#Originally written October 2022, last updated 26/06/2023

import csv

from collections import defaultdict

#Part of the pipeline for DIGITtally analysis - see www.digittally.org
#Analyses expression of orthologs to D. melanogaster genes of interest in the silkworm Bombyx mori using RNAseq data from SilkDB (https://silkdb.bioinfotoolkits.net/main/species-info/-1)


#Gathers tissues of interest from designated text files
def build_tiss_lists(tissuesfile):
    fulllist = set()

    with open(tissuesfile) as INtissues:
        for line in INtissues:
            tissue = line.strip('\n')
            fulllist.add(tissue)

    return fulllist

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

#Populates the B. mori orthologs for each gene of interest, using the designated ortholog file
def pop_orthos(input):
    ortholist = []
    problemlist = {}
    has_an_ortho = set()

    allgenes = []
    dmeltoorthodict = defaultdict(list)
    orthotodmeldict = {}

    with open(input) as foundorthosfile:
        INfos = csv.reader(foundorthosfile)

        next(INfos)
	
	#Each line corresponds to a single gene
        for line in INfos:
            fbgene = line[0].strip('\n')

            allgenes.append(fbgene)

            if line[1] != 'No orthologs found':
                orthos = line[1].strip('[] ')
                if ';' in orthos:
                    indivgenes = orthos.split('; ')
                else:
                    indivgenes = [orthos,]
                
                has_an_ortho.add(fbgene)

                for gene in indivgenes:
                    indivortho = gene.strip(" '")
		    
		    #Trims orthologue names as necessary
                    if indivortho[0].isalpha() == False:
                        indivortho = indivortho[1:]
                    
                    if indivortho[-1].isnumeric() == False:
                        indivortho = indivortho[:-1]
                        
                    #we set up a dictionary of lists which holds all the orthologs associated with each FbID
                    dmeltoorthodict[fbgene].append(indivortho)
                    
                    #The reverse is slightly harder - we want to create a dictionary matching each ortholog to a specific FbID but some orthologs can be associated with multiple
                    if indivortho not in ortholist:
                        orthotodmeldict[indivortho] = fbgene
                        ortholist.append(indivortho)

	            #B.mor genes which are already associated with a FbID are added to the "problem ortho" dictionary
                    else:
                        try:
                            problemlist[indivortho].append(fbgene)
                        except:
                            problemlist[indivortho] = [orthotodmeldict[indivortho], fbgene]
    
    return ortholist, orthotodmeldict, dmeltoorthodict, problemlist, has_an_ortho, allgenes

#Gets the positions of each tissue within each line of the input data
#This is determined at runtime because the SilkDB data is organised in a very unintuitive fashion
def find_tissue_positions(parsed_line):
    itemcount = 0

    tissuetups = defaultdict(list)

    for item in parsed_line:

        if '\n' in item:
            item = item.strip('\n')
            
        if itemcount == 0:
            itemcount += 1

        else:
            subsections = item.split('--')

            tissue = subsections[0]
           
            if '-' in tissue:
                tissue = tissue.replace('-', ' ')
            
            if tissue == 'Fatbody':
                tissue = 'Fat Body'
           
            worm = subsections[1]
            
            indivtuple = (itemcount, tissue)
            tissuetups[worm].append(indivtuple)
            itemcount += 1
    
    return tissuetups

#Creates a list of tuples, of format (expression dictionary, stage name) for each Bombyx mori life stage in use
def create_type_tuples(worm_stages):
    type_tuples = []

    for stage_of_interest in worm_stages:
        indiv_tuple = ({}, stage_of_interest)
        type_tuples.append(indiv_tuple)
    
    return type_tuples

#The tissues which exist in each life stage are added to the appropriate expression dictionary, with the expresssion data accomanying
def populate_individual_dict(type_tuples, tissuetups, parsed_line):

    for individual_type in type_tuples:
        itype_dict = individual_type[0]
        itype_annotation = individual_type[1]

        for item in tissuetups[itype_annotation]:
            itype_dict[item[1]] = parsed_line[item[0]] 

    return type_tuples

#Not all tissues in SilkDB are available for all life stages. 
#Using user-selected tissue types, we check which life stages contain usable information. Only this subset of data will be analysed
def find_usability(types_tuples, targettissues, stages):
    tissues_usable_by_stage = defaultdict(list)
    wormtypes_with_tissue = defaultdict(list)

    refined_targets = [] + list(stages)

    for wormtype in types_tuples:
        
        type_name = wormtype[1]

        for tissue in targettissues:
            try:
                testing = wormtype[0][tissue]
            
                tissues_usable_by_stage[type_name].append(tissue)
                wormtypes_with_tissue[tissue].append(type_name)

            except:
                pass
        
        if tissues_usable_by_stage[type_name] == []:
            del tissues_usable_by_stage[type_name]

            if type_name in refined_targets:
                print(f'None of the target tissues can be found in {type_name} data, so this life stage will not be considered')
                refined_targets.remove(type_name)

    for tissue in wormtypes_with_tissue:
        if wormtypes_with_tissue[tissue] == []:
            print(f'FATAL ERROR - No selected B.mori life stage has data associated with {tissue} tissue. Please run again with default weights')
            exit()

    return tissues_usable_by_stage, wormtypes_with_tissue, refined_targets

#Sets up the following system for marking off individual tissues of interest in which orthologs of a given gene have been found:
    #{Gene:
    #   {Life stage:       
    #       {Category:
    #           [Individual, Tissues]
    #       }
    #   }
    # }
#By removing tissues one at a time when any ortholog of a given Drosophila gene fulfils a criteria, that "category" dictionary will be empty if necessary conditions are fulfilled in all tissues of interest
#Thus, we can track, across various orthologs, which Drosophila genes have orthologs Enriched/Abundant/Expressed in a set of tissues of interest
def initialise_scoring_system(has_an_ortho, tissue_list, stages):
    uncompleted_cats = {}

    for gene in has_an_ortho:
        uncompleted_cats[gene] = {}

        for stage in stages:
            uncompleted_cats[gene][stage] = {}

            for category in ['Enrichment', 'Specificity', 'Expression']:
                uncompleted_cats[gene][stage][category] = [] + list(tissue_list)
    
    return uncompleted_cats

#Creates a dictionary containing details of tissues to be removed from the uncompleted_cats scoring system
#{type:
#   {category:
#       [Individual, Tissues]
#   }
#}
def removal_by_ortholog(types):
    to_remove = {}

    for type in types:
        to_remove[type] = {}

        for category in ['Enrichment', 'Specificity', 'Expression']:
            to_remove[type][category] = []
    
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

#Analyses expression data for orthologs of interest within SilkDB data
def find_silkdb_expression(targettissues, orthologs, orthotodmel, threshold, specthresh, problemorthos, worm_stages, silkdb, has_an_ortho, fpkm_min):
    
    #Creates a list of (expression data dict, silkworm stage) tuples
    base_type_tuples = create_type_tuples(worm_stages)
    
    usability_check = False
    
    with open(silkdb) as silkdbfile:
        linecount = 0

        for line in silkdbfile:
            parsed_line = line.split('\t')
            indivortho = parsed_line[0]
	    
	    #The first line of the data file can be used to find tissue positions 
            if linecount == 0:

                tissuetups = find_tissue_positions(parsed_line)

                linecount += 1
                continue
                	
	    #Each other line in the data contains data on a single gene. 
	    #For the first line , we initialise the scoring system for usable tissue/lifestages
            elif linecount == 1:
                
                #We gather expression data for the gene in all tissue types
                orthologs_types_tuples = populate_individual_dict(base_type_tuples, tissuetups, parsed_line)
                
                #We check which of the tissue types are useable within each worm life stage
                tissues_usable_by_stage, wormtypes_with_tissue, refined_type_list = find_usability(orthologs_types_tuples, targettissues, worm_stages)
		
		#We intialise the uncompleted_cats system for all usable life stages/tissue types
                uncompleted_cats = initialise_scoring_system(has_an_ortho, targettissues, refined_type_list)
                
                linecount += 1
            
            #Other genes are processed similarly
            if indivortho in orthologs:
                abundance_dict = {}

                to_remove = {}

                fbgene = orthotodmel[indivortho]
                
                #We gather expression data for the gene in all tissue types
                orthologs_types_tuples = populate_individual_dict(base_type_tuples, tissuetups, parsed_line)
		
		#Ensures scoring is initialised
                if usability_check == False:

                    tissues_usable_by_stage, wormtypes_with_tissue, refined_type_list = find_usability(orthologs_types_tuples, targettissues, worm_stages)

                    uncompleted_cats = initialise_scoring_system(has_an_ortho, targettissues, refined_type_list)

                    usability_check = True
                    
                #Creates a list in "to_remove" for each life stageof interest
                for type in refined_type_list:
                    to_remove[type] = defaultdict(list)
		
                for indiv_type in orthologs_types_tuples:
                    
                    type_label = indiv_type[1]
		    
		    #For each USABLE worm type of interest:
                    if type_label in refined_type_list:
                        wholeworm = []

                        nontarget_reads = []
			
			#We gather INDIVIDUAL TISSUE expression data
                        reads_in_type = indiv_type[0]

                        for indiv_tissue in reads_in_type:

                            abundance = float(reads_in_type[indiv_tissue])
			    
			    #SilkDB has no true "Whole Insect" - A pseudo-Whole is created by averaging across all tissue types
                            wholeworm.append(abundance)

                            abundance_dict[indiv_tissue] = abundance
                           
                            #We also add the read to a list of non tatget reads if not a target tissue
                            if indiv_tissue not in targettissues:
                                nontarget_reads.append(abundance)
                        
                        #The pseudo-Whole is created by averaging across all tissue types
                        wholewormavg = float(sum(wholeworm)/len(wholeworm))
                        
                        #The pseudo-Whole is adjusted up to the FPKM baseline
                        if wholewormavg < fpkm_min:
                            wholewormavg = fpkm_min

                        for indivtissue in abundance_dict:
                            
                            #For each target tissue, we check for Enrichment, Specificity, and Expression
                            if indivtissue in targettissues and type_label in refined_type_list:

                                enrichment = float(abundance_dict[indivtissue]/wholewormavg)
				
				#If a gene is higher than the pseudo-Whole expression, it is Enriched in that tissue
                                if enrichment > threshold:
                                    to_remove[type_label]['Enrichment'].append(indivtissue)
                                
                                #If a gene is higher than the FPKM threshold, it is Expressed in that tissue
                                if abundance_dict[indivtissue] > fpkm_min:
                                    to_remove[type_label]['Expression'].append(indivtissue)
                                
                                #If a gene is higher than in all non-target tissues, it is Specific to that tissue
                                if abundance_dict[indivtissue] > (max(nontarget_reads) * specthresh):
                                    to_remove[type_label]['Specificity'].append(indivtissue)
		
		#If an multiple FBgns correspond to an ortholog, ALL OF THOSE datasets should mirror.
                if indivortho in problemorthos:
                    for othergene in problemorthos[indivortho]:
                        uncompleted_cats[othergene] = effect_removal(to_remove, uncompleted_cats[othergene])
                
                else:
                    uncompleted_cats[fbgene] = effect_removal(to_remove, uncompleted_cats[fbgene])
    try:
        return uncompleted_cats, wormtypes_with_tissue, tissues_usable_by_stage
    except:
        return {}, wormtypes_with_tissue, tissues_usable_by_stage

#The uncompleted_cats dictionary is resolved by counting empty lists
#If the number of empty lists for each gene, for each category, is equal to the number of sexes of interest, that category scores a "1" in the MozAtlas tally
#Otherwise, the category scores 0
def score_by_cat(uncompleted_cats, tissues_usable_by_stage, target_tissues):

    pos_by_cat = {
        'Enrichment': 0, 'Specificity': 1, 'Expression': 2
    }

    usable_types = [wormtype for wormtype in tissues_usable_by_stage]
    
    usable_types.append('ALL SEARCHABLE TISSUES IN ANY DESIRED TYPE')

    for_all_any = []

    for stage in tissues_usable_by_stage:
        
        if len(tissues_usable_by_stage[stage]) == len(target_tissues):
            for_all_any.append(stage)
    
    if len(for_all_any) > 0:
        usable_types.append('ALL TISSUES OF INTEREST IN ANY DESIRED TYPE')
    
    else:
        print('None of the desired target silkworm life stages have data for ALL the desired target tissues')
        print('Thus, we can only check whether any silkworm type expresses orthologs of a given flybase gene in all the desired tissues it possesses - ie ALL SEARCHABLE TISSUES IN ANY DESIRED TYPE')
        print('Do note that this tends to be a little more promiscuous than the standard check for ALL TISSUES OF INTEREST IN ANY DESIRED TYPE')

    usable_types.append('ALL SEARCHABLE IN ALL DESIRED TYPES')

    tally_compilation = {wormtype : defaultdict(silkdb_tally_default) for wormtype in usable_types}

    for gene in uncompleted_cats:
        
        gene_cats = {'Enrichment': 0, 'Specificity': 0, 'Expression': 0}
        num_types = 0

        for wormtype in uncompleted_cats[gene]:
            num_types += 1

            for category in uncompleted_cats[gene][wormtype]:

                if uncompleted_cats[gene][wormtype][category] == []:
                    gene_cats[category] += 1

                    tally_compilation[wormtype][gene][pos_by_cat[category]] = 1

                    tally_compilation['ALL SEARCHABLE TISSUES IN ANY DESIRED TYPE'][gene][pos_by_cat[category]] = 1

                    if wormtype in for_all_any:
                        tally_compilation['ALL TISSUES OF INTEREST IN ANY DESIRED TYPE'][gene][pos_by_cat[category]] = 1

        for cat in pos_by_cat:

            if gene_cats[cat] == num_types:

                tally_compilation['ALL SEARCHABLE IN ALL DESIRED TYPES'][gene][pos_by_cat[cat]] = 1

    return tally_compilation

#Prints all relevant data to an output file for easy analysis. 
#An accompanying readme file detailing tissue types by ife stage is also generated
def printoutscores(tally_comp, dmeltoortho, threshold, specthresh, outfolder, wormtypes_by_tissue, tissues_usable_by_stage, allgenes):

    with open(f'{outfolder}/Bmor_SilkDB_README.txt', 'w+') as readmefile:
        
        readmefile.write(f'Tissues usable by Bombyx mori life stage\n')

        for wormtype in tissues_usable_by_stage:
            readmefile.write(f'\n{wormtype}:\t{tissues_usable_by_stage[wormtype]}')
        
        readmefile.write(f'\n\nBombyx mori life stage in which a given tissue of interest is found\n')

        for tissue in wormtypes_by_tissue:
            readmefile.write(f'\n{tissue}:\t{wormtypes_by_tissue[tissue]}')

    with open(f'{outfolder}/Bmor_SilkDB_Breakdown.csv', 'w+') as conservefile:
        OUTconserve = csv.writer(conservefile)
        
        #This block builds the header row
        header = ['FB_ID', 'Bombyx mori ortholog(s)']

        for wormtype in tally_comp:
            if len(wormtype) == 1:
                header_extension = [f'Any Ortholog combination enriched >= {threshold} in target tissues vs whole mosquito - {wormtype} insects', f'Any Ortholog combination {specthresh} times higher abundance in target tissues compared to non target tissues? - {wormtype} insects', f'ANY orthologs expressed in all target tissues? - {wormtype} insects']
                header.extend(header_extension)

            else:
                header_extension = [f'{wormtype} - enrichment > {threshold}', f'{wormtype} - increased abundance over non-targets', f'{wormtype} - any expression']
                header.extend(header_extension)
        
        OUTconserve.writerow(header)

        #Gene data is printed: [D.mel gene, A.gam orthos (if any), then enrichment > threshold, abundance > non-target tissues, expression  in tissues of interest PER FLY TYPE]
        for gene in allgenes:

            try:
                if len(dmeltoortho[gene]) > 1:
                    bmorthologs = '; '.join(dmeltoortho[gene])
                else:
                    bmorthologs = dmeltoortho[gene][0]
            except:
                bmorthologs = 'No orthologs found'

            outrow = [gene, bmorthologs, ]

            for wormtype in tally_comp:
                outrow.extend(tally_comp[wormtype][gene])

            OUTconserve.writerow(outrow)

#Creates the default tally for SilkDB
def silkdb_tally_default():
    return [0, 0, 0]

#Creates a list of silkworm types which the user wishes included
def generate_use_list(weight_dict):
    use_list = []

    for category in weight_dict:
        if weight_dict[category] > 0:
            use_list.append(category)
    
    return use_list

#This function uses the user-supplied weighting for the various readings (Enrichment, Abundance, Expression) and for MozTubules (compared to a tissue from MozAtlas)
#The final tally will be a score OUT OF ONE based on the entirity of Anopheles gambiae orthology
def bmor_to_digittally(allgenes, tally_compilation, measure_weights, lifestage_weights, weight_dict, types_in_use):

    #This will be a list of lists for csv printing
    gene_rows = []

    #The divisor_by_type generates a per-fly-type possible maximum score for MozAtlas
    divisor_by_type = (measure_weights['Enrichment'] + measure_weights['Specificity'] +  measure_weights['Expression'])

    categories_to_check = []

    individual_weights = {}

    if lifestage_weights['Individual_stage'] != 0:
        for stage in weight_dict:
            if weight_dict[stage] != 0 and stage in types_in_use:
                categories_to_check.append(stage)
                individual_weights[stage] = weight_dict[stage]
    
    if lifestage_weights['All_type'] != 0:
        categories_to_check.append('ALL SEARCHABLE IN ALL DESIRED TYPES')
        individual_weights['ALL SEARCHABLE IN ALL DESIRED TYPES'] = lifestage_weights['All_type']
    
    if lifestage_weights['Any_type'] != 0:

        try:
            test = tally_compilation['ALL TISSUES OF INTEREST IN ANY DESIRED TYPE']
            categories_to_check.append('ALL TISSUES OF INTEREST IN ANY DESIRED TYPE')
            individual_weights['ALL TISSUES OF INTEREST IN ANY DESIRED TYPE'] = lifestage_weights['Any_type']

        except:
            categories_to_check.append('ALL SEARCHABLE TISSUES IN ANY DESIRED TYPE')
            individual_weights['ALL SEARCHABLE TISSUES IN ANY DESIRED TYPE'] = lifestage_weights['Any_type']


    for gene in allgenes:

        total_divisor = 0
        silkdb_score = 0

        #For each fly type, a score is compiled based on the weights of each ranking component.
        #Effectively, this score = the weight for each component if the requirements are met.
        #The total divisor = the divisor_by_type times the number of wormtypes
        for wormtype in categories_to_check:
            silkdb_score += float(measure_weights['Enrichment'] * tally_compilation[wormtype][gene][0]) * individual_weights[wormtype]
            silkdb_score += float(measure_weights['Specificity'] * tally_compilation[wormtype][gene][1]) * individual_weights[wormtype]
            silkdb_score += float(measure_weights['Expression'] * tally_compilation[wormtype][gene][2]) * individual_weights[wormtype]

            total_divisor += divisor_by_type * individual_weights[wormtype]
        
        #The final silkdb score = the total score across desired flytypes over the total divisor - ie, reduced to a 0-1 scale
        final_silkdb_score = silkdb_score/total_divisor

        #A row to be printed for each gene is then stored in gene_rows
        gene_rows.append([gene, final_silkdb_score])
    
    return gene_rows

#This function creates the "Agam_orthology_tally csv, needed for the final DIGITtally output"
def append_to_tally(gene_rows, outfolder):
    with open(f'{outfolder}/Bmor_orthology_tally.csv', 'w+') as tallyfile:

        OUTtally = csv.writer(tallyfile)

        header = ['Gene symbol', 'Bombyx mori score']
        OUTtally.writerow(header)

        for line in gene_rows:
            OUTtally.writerow(line)

def Score_Bmor_orthology(ortholog_file, enrichthresh, specthresh, tissues, silkdb_file, outputorthology, outputtally, fpkm_background, typeweights, lifestageweight, measureweights):

    target_tissues = tissues

    #We read the user input weights into a dictionary (Key = Fly type, value = weight/total weights)
    weight_dict = gather_weights(typeweights)
    measure_weights = gather_weights(measureweights)
    lifestage_weights = gather_weights(lifestageweight)
    worm_stages_in_use = generate_use_list(weight_dict)

    orthologs, orthotodmel, dmeltoortho, problematicorthos, has_an_ortho, all_genes = pop_orthos(ortholog_file)

    categories, wormtypes_by_tissue, tissues_usable_by_stage = find_silkdb_expression(target_tissues, orthologs, orthotodmel, enrichthresh, specthresh, problematicorthos, worm_stages_in_use, silkdb_file, has_an_ortho, fpkm_background)

    tally_compilation = score_by_cat(categories, tissues_usable_by_stage, target_tissues)

    printoutscores(tally_compilation, dmeltoortho, enrichthresh, specthresh, outputorthology, wormtypes_by_tissue, tissues_usable_by_stage, all_genes)

    #The DIGITtally output is ALWAYS produced
    rows_to_write = bmor_to_digittally(all_genes, tally_compilation, measure_weights, lifestage_weights, weight_dict, tissues_usable_by_stage)
    append_to_tally(rows_to_write, outputtally)
