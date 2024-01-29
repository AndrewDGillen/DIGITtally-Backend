import csv

from collections import defaultdict

#Creates a list of the desired annotations, from the annotation text file.
def populate_lists(ids):

    idlist = set()

    with open(ids) as idfile:
        for line in idfile:
            idlist.add(line.strip('\n'))
    
    return idlist

#Uses the FlyBase Human orthology data to populate the orthologs associated with each gene of interest
#All hits with a DIOPT score over 3 by default (Moderate or better) are listed as possible orthologs for a biven FB gene
#The best match for each gene, the associated diseases, and score is also found
def find_orthologs(genelist, fborthos, stringency):
    databygene = {}

    #A defaultdict is created for each category of data, which will contain a set (for ease of searching)
    for cat in ['PluralNames', 'PluralDiseases', 'MonoName', 'MonoDiseases', 'MonoScore', 'Scores']:
        databygene[cat] = defaultdict(set)

    with open(fborthos) as orthofile:

        INdioptdata = csv.reader(orthofile, delimiter = '\t')

        for line in INdioptdata:

            #These IF statements prevent processing of inappropriate lines
            #This is structured in a somewhat ugly way to futureproof against changes in the number of annotation lines in the DIOPT file
            if line != []:
                if line[0][0] != '#':

                    #Data associated with each ortholog is harvested
                    #Ortholog relationships are ONE PER LINE - ie each flybase gene may have information across multiple lines
                    fbid = line[0]
                    orthoname = line[4]
                    score = int(line[5])
                    diseases = line[7]

                    #If the gene is a gene of interest, we store the relevant info in our data dictionary
                    if fbid in genelist:
                        
                        if score >= stringency:
                            databygene['PluralNames'][fbid].add(orthoname)
                            databygene['PluralDiseases'][fbid].add((diseases))
                            databygene['Scores'][fbid].add((score, orthoname, diseases))
        
        #We then sort the data for each gene to: 
        #   -Check whether each gene has an ortholog associated
        #   -Find the best ortholog match for each gene
        for gene in genelist:
            
            #If a given gene has no associated score, it has no known human orthologs, so an empty set is returned
            if gene not in databygene['Scores']:
                for cat in ['PluralNames', 'PluralDiseases', 'Scores']:
                    databygene[cat][gene] = {''}

            #We produce a list of scores to be sorted to find the highest scoring ortholog for each gene
            sorteddata = list(databygene['Scores'][gene])
            
            #This carried out the sort
            try:
                sorteddata.sort(key=lambda i:i[0],reverse=True)
                databygene['MonoScore'][gene] = sorteddata[0][0]
                databygene['MonoName'][gene] = sorteddata[0][1]
                databygene['MonoDiseases'][gene] = sorteddata[0][2]
            
            #If the sort fails due to being carried out on an empty list, "N/A" is returned
            except:
                databygene['MonoScore'][gene] = 'N/A'
                databygene['MonoName'][gene] = 'N/A'
                databygene['MonoDiseases'][gene] = 'N/A'
    
    return databygene

#A file containing a full summary of Human Orthology for each D. mel gene of interest
def make_orthofile(genelist, databygene, outfolder, stringency):

    with open(f'{outfolder}/Human_Ortholog_DIOPTResults.csv', 'w+') as orthofile:
        OUTorthos = csv.writer(orthofile)

        header = ['FlyBase ID', f'All orthologs with a DIOPT score > {stringency}', 'Human Diseases associated with these genes (OMIM phenotype ID[Full name])', ' ', 'Highest scoring ortholog', 'Highest Score', 'Disease associated with the best ortholog match']
        OUTorthos.writerow(header)
        
        #An output row is created for each gene
        for gene in genelist:
            outrow = [gene,]

            for cat in ['PluralNames', 'PluralDiseases', ' ', 'MonoName', 'MonoScore', 'MonoDiseases']:

                #This just maintains correct output spacing
                if cat == ' ':
                    outrow.append(cat)
                
                #The categories with multiple entries are handled in a specific manner to make lists more readable
                elif cat in ['PluralNames', 'PluralDiseases']:

                    #Try/Except catches genes with no databygene entry
                    try:
                        #We convert the lists to a more legible string for printing
                        printable = ', '.join(databygene[cat][gene])
                        
                        #If this string is empty, we instead report that no ortholgs/diseases have been found
                        if printable == ', ' or printable == '':
                            printable = 'None found'
                        
                        #This strips the often strangely-formatted lists of associated diseases down to an easily legible format
                        else:
                            while printable[0] == ',' or printable[-1] == ' ':
                                printable = printable.strip(', ')

                        outrow.append(printable)

                    except:
                        outrow.append('None found')
                
                #The mono entry categories can be printed as is
                else:
                    outrow.append(databygene[cat][gene])

            OUTorthos.writerow(outrow)

#This function creates the "Hsap_orthology_tally csv, needed for the final DIGITtally output"
def make_orthotally(genelist, databygene, outfolder, detected_weight, disease_weight):

    with open(f'{outfolder}/Hsap_orthology_tally.csv', 'w+') as tallyfile:
        OUTtally = csv.writer(tallyfile)

        header = ['FlyBase ID', 'Human Orthology Score']

        OUTtally.writerow(header)
        
        for gene in genelist:
            outrow = [gene,]
            
            tally = 0
            target = max(detected_weight, 0) + max(disease_weight, 0)

            #If the program is run, we always want to include whether a gene has a human ortholog in the H. sapiens score.
            if databygene['PluralNames'][gene] != []:
                tally += (1*detected_weight)

            #If we're interested in whether a human ortholog is associated with a human disease, we add to the tally if the "plural diseases" category is not empty.
            #Then, since we want the Hsap score to be out of one, we divide the total tally by 2.
            if disease_weight != 0:
                if databygene['PluralDiseases'][gene] != {''}:
                    tally += (1*disease_weight)

            tally = tally/target

            outrow.append(tally)

            OUTtally.writerow(outrow)

def AnalyseDIOPTFindings(input_fbids, detection_weight, disease_weight, stringency, outputorthos, outputtally, flybase_file):

    goilist = populate_lists(input_fbids)
    databygene = find_orthologs(goilist, flybase_file, stringency)

    make_orthofile(goilist, databygene, outputorthos, stringency)
    make_orthotally(goilist, databygene, outputtally, detection_weight, disease_weight)
