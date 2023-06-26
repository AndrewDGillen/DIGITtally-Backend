#Code written by Dr. Andrew Gillen (Dow-Davies lab, University of Glasgow)
#Originally written October 2022, last updated 26/06/2023

import csv
import gc
import loompy as lp
import pandas as pd

#Part of the pipeline for DIGITtally analysis - see www.digittally.org
#Specifically, this program checks whether genes of interest, or any known synonym, can be identified within FlyCellAtlas (https://flycellatlas.org/)

#Creates lists of genes of interest and reference genes
def populate_list(gois):

    goilist = []

    with open(gois) as goifile:
        for line in goifile:
            goilist.append(line.strip('\n'))
    
    return goilist

#takes incoming FlyBase IDs and converts them to current Symbols for SCOPE.
def convert_refgenes(references, synonymfile):
    formatted_reflist = []

    with open(synonymfile) as synfile:  
        INsyns = csv.reader(synfile, delimiter = '\t')

        for i in range(6):
            next(INsyns)

        for line in INsyns: 
            if line[0] in references and line[0] != '':
                print(line[0],line[2])
                formatted_reflist.append(line[2])

    print(formatted_reflist)
    return(formatted_reflist)

#Uses the FlyCellAtlas LOOM file to check that genes in the supplied list are recognised by the FCA dataset
def firstpass(genelist, df):
    genesfound = []
    notfound = []

    for indivgene in genelist:
        try:
            #We don't actually do anything with this dataframe, but trying to create it will raise an exception if the gene isn't in the loom file. 
            #This is the easiest way I've found to catch genes not in FCA.
            geneframe = pd.DataFrame(df[df.ra.Gene == indivgene,:], index=[indivgene], columns=df.ca["CellID"]).T

            genesfound.append(indivgene)

            #Frees up some memory as these are big frames
            del geneframe
            gc.collect()
            geneframe = pd.DataFrame()

        #Genes which can't be used to form the dataframe are returned in the notfound list
        except:
            notfound.append(indivgene)

    return genesfound, notfound

#Uses the FlyBase_Synonyms file to find known synonyms for genes which cannot be found in SCOPE by default
def find_synonyms(notfound, synonymfile):
    synonymdict = {}
    unfound = []

    with open(synonymfile) as synfile:  
        INsyns = csv.reader(synfile, delimiter = '\t')

        for line in INsyns:
            try:
                #Synonyms are found | seperated at position 5 in each line
                synonyms = line[5].split('|')

                if line != '':
                    #FlyBase_Synonyms is a BIG file containing synonyms for many different IDs in different species
                    #For our purposes, we only consider genes (FBgn) from Drosophila melanogaster (Dmel)
                    if line[0][:4] == 'FBgn' and line[1] == 'Dmel':

                        for item in notfound:
                            
                            #If the known symbol is found, we assign each of its synonyms to the synonymdict
                            if item == line[2]:
                                
                                if item not in synonymdict:
                                    synonymdict[item] = []
                                
                                synonymdict[item] = synonymdict[item] + synonyms
            except:
                pass
    
    #Genes for which synonyms CANNOT be found are appended to the "unfound" list
    for gene in notfound:
        if gene not in synonymdict:
            unfound.append(gene)
    
    return synonymdict, unfound

#Rechecks the synonyms in the same manner as before. 
def recheck_synonyms(synonymdict, unfindable, genesfound, df):
    foundsynonyms = []

    for gene, synonyms in synonymdict.items():
        foundcheck = 0
        
        #Each individual synonym for a given gene is tried before giving up
        for indivsyn in synonyms:
            
            #strips off newlines before we check the dataframe as these aren't tolerated
            indivsyn = indivsyn.strip('\n')

            try:
                geneframe = pd.DataFrame(df[df.ra.Gene == indivsyn,:], index=[indivsyn], columns=df.ca["CellID"]).T

                genesfound.append(indivsyn)
                foundsynonyms.append([gene, indivsyn])
                foundcheck += 1
                break

            except:
                pass
        
        #In this case, if NO synonym for a gene can be found in FCA, we have to accept that we cannot gather data for that gene from the FCA dataset easily.
        if foundcheck == 0:
            unfindable.append(gene)
    
    return foundsynonyms, genesfound, unfindable

#Creates "Found checked" and "unfound" gene lists, simple .txt files
def make_list_output(desiredname, suppliedlist):
    with open(desiredname, 'w+') as OUTlist:
        for item in suppliedlist:
            OUTlist.write(f'{item}\n')

#Alters "synonym" .csv, created earlier in the pipeline, which can be used to match up original gene names with the synonym used for the FCA search.
def make_synonym_output(desiredname, suppliedlist):
    with open(desiredname, 'w+') as synonymfile:
        OUTsynonym = csv.writer(synonymfile)

        header = ['Symbol in DetectedGOIs_Symbol.txt', 'UPDATED Symbol for DetectedGOIs_Symbol_Checked.txt']
        OUTsynonym.writerow(header)

        for synpair in suppliedlist:
            OUTsynonym.writerow(synpair)

#A big function whihc runs all the necessary components and keeps the end user up to speed.
def full_check_recheck(genelist, genetype, outfolder, df, synonymfile):
    goitarget = len(genelist)

    genesfound, notfound = firstpass(genelist, df)

    #Only rechecks if any genes are actually missed
    if len(genesfound) < goitarget:
        print(f'\nNot all genes immediately recognized by SCOPE. {len(notfound)} will be checked for known synonyms.')

        synonymdict, unfindable = find_synonyms(notfound, synonymfile)
        
        oldunfindlen = len(unfindable)

        if len(unfindable) > 0:
            print(f'\nNot all genes are recognized by FlyBase. {len(unfindable)} will be excluded from FlyCellAtlas analysis')

        foundsynonyms, genesfound, unfindable = recheck_synonyms(synonymdict, unfindable, genesfound, df)

        if len(unfindable) > oldunfindlen:
            print(f'\nNot all of the genes have a synonym recognised by SCOPE. In total, {len(unfindable)} will be excluded from FlyCellAtlas analysis.\n Please see the NotInSCOPE genelist for more information.')
        
        make_list_output(f'{outfolder}/{genetype}_Symbol_checked.txt', genesfound)
        make_list_output(f'{outfolder}/{genetype}_Symbol_NotInSCOPE.txt', unfindable)
        make_synonym_output(f'{outfolder}/Synonyms/{genetype}_Synonyms.csv', foundsynonyms)

        print(f'\nIn total, {len(genesfound)} {genetype} (of {len(genelist)} checked) have a symbol which can be detected in SCOPE. These can be used to proceed to FlyCellAtlas analysis')

    else:
        make_list_output(f'{outfolder}/{genetype}_Symbol_checked.txt', genesfound)
        print(f'\nAll desired {genetype} found on the first pass. These can be used to proceed to FlyCellAtlas analysis')

def CheckTheList(loomfile, genes_of_interest, references, output_folder, synonymfile):

    goilist = populate_list(genes_of_interest)
    
    df = lp.connect(loomfile, validate = False)

    full_check_recheck(goilist, 'DetectedGOIs', output_folder, df, synonymfile)

    if len(references) > 0:
        reflist = convert_refgenes(references, synonymfile)
        full_check_recheck(reflist, 'UserRefGenes', output_folder, df, synonymfile)
