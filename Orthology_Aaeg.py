#Code written by Dr. Andrew Gillen (Dow-Davies lab, University of Glasgow)
#Originally written October 2022, last updated 26/06/2023

import csv

from collections import defaultdict

#Part of the pipeline for DIGITtally analysis - see www.digittally.org
#Analyses expression of orthologs to D. melanogaster genes of interest in the mosquito Aedes Aegypti using RNAseq data from Aegypti-Atlas (http://aegyptiatlas.buchonlab.com/)

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

#Populates the A. aeg orthologs for each gene of interest, using the designated ortholog file
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
	
	#Each line corresponds to a single gene
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

                    #A.aeg genes which are already associated with a FbID are added to the "problem ortho" dictionary
                    else:
                        
                        #This ensures correct layout of the dict - {A.aeg gene: [list of flybase IDs]}
                        try:
                            problem_orthos[gene].append(fbgene)
                        except:
                            problem_orthos[gene] = [orthotodmeldict[gene], fbgene]
    
    return ortholist, orthotodmeldict, dmeltoorthodict, problem_orthos, has_an_ortho, all_genes

#Gets the average expression of a gene in a single tissue, given the indices which correspond to that tissue
def get_tissue_expression(indices, current_line):
    expression = float(0)

    for ind_index in indices:
        expression += float(current_line[ind_index])
    
    expression = float(expression / len(indices))
    
    return expression

#Just a clean, readable layout to produce the tissue indices dictionary, containing the positions of each read associated with a given tissue in AegyptiAtlas
#Tissues are defined either as "Basic" or as "Gut Subset" to aid in defining non-target tissues
def define_tissues():
    tissue_indices = {
        'Head' : ([10, 12, 14], 'Basic'),
        'Thorax' : ([16, 18, 20], 'Basic'),
        'Abdomen' : ([22, 24, 26], 'Basic'),
        'Gut (Generic)' : ([28, 30, 32, 34, 36, 38], 'Basic'),
        'Malpighian Tubules' : ([40, 42, 44], 'Basic'),
        'Ovaries' : ([46, 48, 50], 'Basic'),
        'Crop' : ([54, 56, 58, 60], 'Gut subset'),
        'Proventriculus' : ([62, 64, 66], 'Gut subset'),
        'Anterior Midgut' : ([68, 70, 72], 'Gut subset'),
        'Posterior Midgut' : ([74, 76, 78], 'Gut subset'),
        'Hindgut' : ([80, 82, 84], 'Gut subset')
    }

    return tissue_indices

#Uses the user-defined tissues to prepare a list of non-target tissues
#If the user is interested in a gut subset, other gut subsets will be added to this list, but not the generic form of Gut
#If not, the gut subsets will not be added
#This prevents, for example, genes enriched in the Generic Gut tissue from being excluded if enriched in a specific subset of gut tissues.
def define_non_targets(tissue_list, tissue_indices):
    
    use_gut_subset = False
    tissues_for_nontarget = []

    for indiv_tissue in tissue_list:
        print('desired:', indiv_tissue)
        if tissue_indices[indiv_tissue][1] == 'Gut subset':
            use_gut_subset = True
    
    if use_gut_subset == True:

        for tissue in tissue_indices:

            if tissue != 'Gut (generic)':
                tissues_for_nontarget.append(tissue)
    
    else:

        for tissue in tissue_indices:
            
            if tissue_indices[tissue][1] == 'Basic':
                tissues_for_nontarget.append(tissue)
    print(tissues_for_nontarget)
    return tissues_for_nontarget

#Sets up the following system for marking off individual tissues of interest in which orthologs of a given gene have been found:
    #{Gene:     
    #   {Category:
    #       [Individual, Tissues]
    #   }
    # }
#By removing tissues one at a time when any ortholog of a given Drosophila gene fulfils a criteria, that "category" dictionary will be empty if necessary conditions are fulfilled in all tissues of interest
#Thus, we can track, across various orthologs, many of which have multiple Microarray probes associated, which Drosophila genes have orthologs Enriched/Abundant/Expressed in a set of tissues of interest
def initialise_scoring_system(has_an_ortho, targettiss):
    uncompleted_cats = {}

    for gene in has_an_ortho:
        uncompleted_cats[gene] = {}

        for category in ['Enrichment', 'Specificity', 'Expression']:
            uncompleted_cats[gene][category] = [] + list(targettiss)
    
    return uncompleted_cats

#Uses the "to_remove" dictionary to remove individual tissues for which conditions have been fulfilled from the "Uncompleted_cats" dictionary on a BY GENE basis
#Notably, this can be called on genes which have already been processed previously, removing tissues in which conditions may not have been met on the previous pass
#THIS IS A VITAL FEATURE
def effect_removal(to_remove, uncompleted_cats_section):
    newcats = uncompleted_cats_section
    
    for category in to_remove:

        for tissue in to_remove[category]:

            if tissue in newcats[category]:

                newcats[category].remove(tissue)
    
    return newcats

#Gets the expression of a list of genes from AegyptiAtlas data
def find_aegyptiatlas_expression(tissues_of_interest, orthologs, orthotodmel, threshold, specthresh, atlasdata, fpkm_background, problemorthos, has_an_ortho):

    #Defines target/Non-target tissues
    tissue_indices = define_tissues()
    tissues_to_check = define_non_targets(tissues_of_interest, tissue_indices)

    #Sets up uncompleted_cats dictionary, which keeps track of which tissues express which gene
    uncompleted_cats = initialise_scoring_system(has_an_ortho, tissues_of_interest)

    with open(atlasdata) as aeatlasfile:
        INaareader = csv.reader(aeatlasfile)
	
	#Skips past headers
        for i in range(3):
            next(INaareader)

	#Each line corresponds to a single gene
        for line in INaareader:
            
            #If the gene is an ortholog of interest, we analyse its expression 
            if line[0] in orthologs:
                enrichment = {}
                expression = {}

                indivortho = line[0]
                fbgene = orthotodmel[indivortho]

                to_remove = defaultdict(list)

                nontarget_expression = []

                #Whole body is handled seperately due to the need to set a lower boundary
                wholebody = get_tissue_expression([4, 6, 8], line)

                #Sets a lower boundary for whole body expression to control enrichment readings
                if wholebody < fpkm_background:
                    wholebody = fpkm_background

                #Expression and enrichment values are gathered for the tissues we want to check
                for tissue in tissues_to_check:

                    current_expression = get_tissue_expression(tissue_indices[tissue][0], line)

                    if tissue not in tissues_of_interest:
                        nontarget_expression.append(current_expression)
                    
                    else:
                        expression[tissue] = current_expression

                        #As the baseline for expression in this dataset is one, if a
                        if current_expression > fpkm_background:
                            to_remove['Expression'].append(tissue)

                        enrichment = expression[tissue] / wholebody

                        if enrichment > threshold:
                            to_remove['Enrichment'].append(tissue)

                for tissue in tissues_of_interest:

                    if expression[tissue] > fpkm_background and expression[tissue] > (max(nontarget_expression) * specthresh):
                        to_remove['Specificity'].append(tissue)

                if indivortho in problemorthos:
                    for othergene in problemorthos[indivortho]:

                        uncompleted_cats[othergene] = effect_removal(to_remove, uncompleted_cats[othergene])
                    
                else:
                    uncompleted_cats[fbgene] = effect_removal(to_remove, uncompleted_cats[fbgene])

    return uncompleted_cats

#The uncompleted_cats dictionary is resolved by counting empty lists
#If the number of empty lists for each gene, for each category, is equal to the number of sexes of interest, that category scores a "1" in the MozAtlas tally
#Otherwise, the category scores 0
def score_by_cat(uncompleted_cats, mozatlas_tally):
    
    pos_by_cat = {
        'Enrichment': 0, 'Specificity': 1, 'Expression': 2
    }
        
    for gene in uncompleted_cats:
        
        gene_cats = {'Enrichment': 0, 'Specificity': 0, 'Expression': 0}

        for category in uncompleted_cats[gene]:

            if uncompleted_cats[gene][category] == []:
                gene_cats[category] += 1

                mozatlas_tally[gene][pos_by_cat[category]] = 1
        
    return mozatlas_tally

#Prints relevant gene expression data to an output file
def printoutscores(hitscores, dmeltoortho, threshold, specthresh, outpath, all_genes):

    with open(f'{outpath}/Aaeg_AegyptiAtlas_Breakdown.csv', 'w+') as conservefile:
        OUTconserve = csv.writer(conservefile)
        
        header = ['Gene Symbol', 'A. Aegypti Ortholog(s)', f'Any combination of orthologs enriched in all target tissues > {threshold}?', 
        f'Any combination of orthologs {specthresh} more abundant in target tissues than non-targets', 'Any combination of orthologs expressed in all target tissues']
                
        OUTconserve.writerow(header)

        count = 0
        for gene in all_genes:
            count += 1
            try:
                if len(dmeltoortho[gene]) > 1:
                    agorthologs = '; '.join(dmeltoortho[gene])
                else:
                    agorthologs = dmeltoortho[gene][0]
            except:
                agorthologs = 'No orthologs found'

            outrow = [gene, agorthologs, ]
            outrow.extend(hitscores[gene])

            OUTconserve.writerow(outrow)

#Creates the default tally for AegyptiAtlas
def aegypti_atlas_default():
    return [0, 0, 0]

#This function uses the user-supplied weighting for the various readings (Enrichment, Abundance, Expression) and for MozTubules (compared to a tissue from MozAtlas)
#The final tally will be a score OUT OF ONE based on the entirity of Anopheles gambiae orthology
def aaeg_to_digittally(allgenes, aegypti_atlas_tally, measure_weights):

    #This will be a list of lists for csv printing
    gene_rows = []

    #The divisor_by_type generates a per-fly-type possible maximum score for MozAtlas
    divisor = (measure_weights['Enrichment'] + measure_weights['Specificity'] +  measure_weights['Expression'])

    for gene in allgenes:
        aegyptiatlas_score = 0

        #For each fly type, a score is compiled based on the weights of each ranking component.
        #Effectively, this score = the weight for each component if the requirements are met.
        #The total divisor = the divisor_by_type times the number of flytypes
        aegyptiatlas_score += float(measure_weights['Enrichment'] * aegypti_atlas_tally[gene][0])
        aegyptiatlas_score += float(measure_weights['Specificity'] * aegypti_atlas_tally[gene][1])
        aegyptiatlas_score += float(measure_weights['Expression'] * aegypti_atlas_tally[gene][2])
            
        #The final mozatlas score = the total score across desired flytypes over the total divisor - ie, reduced to a 0-1 scale
        final_aegyptiatlas_score = aegyptiatlas_score / divisor

        #A row to be printed for each gene is then stored in gene_rows
        gene_rows.append([gene, final_aegyptiatlas_score])
    
    return gene_rows

#This function creates the "Aaeg_orthology_tally csv, needed for the final DIGITtally output"
def append_to_tally(gene_rows, outfolder):
    with open(f'{outfolder}/Aaeg_orthology_tally.csv', 'w+') as tallyfile:

        OUTtally = csv.writer(tallyfile)

        header = ['Gene symbol', 'Anopheles gambiae score']
        OUTtally.writerow(header)

        for line in gene_rows:
            OUTtally.writerow(line)

def Score_Aaeg_Orthology(ortholog_file, enrichthresh, specthresh, tissues, aegyptiatlas_file, outputorthology, outputtally, fpkm_background, measureweights):
    
    target_tissues = tissues

    measure_weights = gather_weights(measureweights)

    orthologs, orthotodmel, dmeltoortho, problematicorthos, has_an_ortho, all_genes = pop_orthos(ortholog_file)

    aegypti_atlas_tally = defaultdict(aegypti_atlas_default)

    categories = find_aegyptiatlas_expression(target_tissues, orthologs, orthotodmel, enrichthresh, specthresh, aegyptiatlas_file, fpkm_background, problematicorthos, has_an_ortho)

    aegypti_atlas_tally = score_by_cat(categories, aegypti_atlas_tally)

    printoutscores(aegypti_atlas_tally, dmeltoortho, enrichthresh, specthresh, outputorthology, all_genes)

    #The DIGITtally output is produced
    rows_to_write = aaeg_to_digittally(all_genes, aegypti_atlas_tally, measure_weights)
    append_to_tally(rows_to_write, outputtally)
