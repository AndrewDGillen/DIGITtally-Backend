#Code written by Dr. Andrew Gillen (Dow-Davies lab, University of Glasgow)
#Originally written October 2022, last updated 26/06/2023

import csv

#Part of the pipeline for DIGITtally analysis - see www.digittally.org
#Gathers Orthologues for DIGITtally species (Aedes aegypti, Anopheles gambiae, Bombyx mori) from OrthoDB data (https://www.orthodb.org/)

#creates ortholog output files for each species and prints the header row. 
def create_file(species, outfolder):

    for_file_name = species.replace(' ', '_')

    with open(f'{outfolder}/{for_file_name}_orthologs.csv', 'w+') as orthologfile:
        OUTorthos = csv.writer(orthologfile)

        header = ['Gene Name', 'Orthologs']
        OUTorthos.writerow(header)

#Simple function for parsing genes of interest
def populate_lists(annos):

    annolist = []

    with open(annos) as annofile:
        for line in annofile:
            annolist.append(line.strip('\n'))
    
    return annolist

#Returns the coluimn of the data file which contains orthologs from a specified species
def get_species_column(species):
    if species == 'Aedes aegypti':
        return 5
    elif species == 'Anopheles gambiae':
        return 4
    elif species == 'Bombyx mori':
        return 8
    else:
        print(f'{species} not recognised in orthology file - please check settings!!')
        exit()

#Returns the orthologs for a list of genes from the data file
def locate_orthologs(species, genes, fullfile):
    
    orthologs = {}
    species_appropriate_column = get_species_column(species)

    with open(fullfile) as all_orthos:
        INorthos = csv.reader(all_orthos)

        next(INorthos)

        for line in INorthos:

            if line[0] in genes:
                iortho = line[species_appropriate_column]

                if iortho == '':
                    orthologs[line[0]] = ['No orthologs found',]
                elif ';' in iortho:
                    orthologs[line[0]] = iortho.split(';')
                else:
                    orthologs[line[0]] = [iortho,]

    return orthologs

#Amends the ortholog files with discovered information
def update_file(species, outfolder, all_gene_orthos):
    for_file_name = species.replace(' ', '_')

    with open (f'{outfolder}/{for_file_name}_orthologs.csv', 'a+') as orthologfile:
        OUTorthos = csv.writer(orthologfile)

        for gene in all_gene_orthos:

            if len(all_gene_orthos[gene]) > 1:
                ortho_to_print = "; ".join(all_gene_orthos[gene])

            elif len(all_gene_orthos[gene]) == 1:
                ortho_to_print = ''.join(all_gene_orthos[gene])

            OUTorthos.writerow([gene, ortho_to_print])

#THIS SERVES AS A REPLACEMENT TO THE API IMPLEMENTATION - Notionally this should speed up DIGITtally to a huge degree
#THE ORHTODB FILE WILL NEED TO BE KEPT UP TO DATE. This meant UPDATES EVERY 2 MONTHS to maintain parity to FlyBase
def gather_orthologs(input_fbids, ortholog_file, species_list, outfolder):

    #The genelist is always populated
    genelist = populate_lists(input_fbids)

    for indiv_species in species_list:
        create_file(indiv_species, outfolder)

        species_orthos = locate_orthologs(indiv_species, genelist, ortholog_file)

        update_file(indiv_species, outfolder, species_orthos)

