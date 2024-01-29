import csv
import glob
import os

import scipy.stats as sp

#Creates lists of tissues/cell types of interest
def populate_list(target, targetlist):

    with open(target) as targetfile:
        for line in targetfile:
            targetlist.append(line.strip('\n'))
    
    return targetlist

# #Creates a dictionary assigning each gene which was changed by SCOPEscope to their original name
# def populate_synonyms(synfile):
#     synonymdict = {}

#     with open(synfile) as synonyms:
#         INsyn = csv.reader(synonyms)
#         next(INsyn)

#         for line in INsyn:
#             synonymdict[line[1]] = line[0]

#Gets the actual designations for target tissues from the FCA_Cell_Types file, along with total counts for target and non-target tissues to be used in statistics.
def get_celltypes(targets, associatedfiles):
    targettypelist= []
    nontargetlist = []

    with open(f'{associatedfiles}/FCA_Cell_Types.csv') as celltypesfile:
        INtypes = csv.reader(celltypesfile)
        next(INtypes)

        tot_dict = {}

        for category in ['Target Tissues', 'ALL Non Target']:
            tot_dict[category] = 0
        
        for line in INtypes:
            acttissname = line[1]
            letsmod = acttissname.split(' : ')
            upper_tiss_type = letsmod[0].capitalize()

            if "_" in upper_tiss_type:
                upper_tiss_type = upper_tiss_type.replace('_', ' ')
            upper_cell_type = letsmod[1].capitalize()
            acttissname = upper_tiss_type + ' : ' + upper_cell_type
            totalcount = line[5]

            if acttissname in targets:
                targettypelist.append(acttissname)
                tot_dict['Target Tissues'] += float(totalcount)

            else:
                nontargetlist.append(acttissname)
                tot_dict['ALL Non Target'] += float(totalcount)

    return targettypelist, nontargetlist, tot_dict

#Frequencies used to check distribution for enrichment within Target Tissues specifically
def get_freqs(tot_dict):
    totcount_t = tot_dict['Target Tissues']
    totcount_nt = tot_dict['ALL Non Target']

    totalcount = totcount_t + totcount_nt

    freq_te = float(totcount_t/totalcount)
    freq_nt = float(totcount_nt/totalcount)

    return freq_te, freq_nt

#desc
def handle_positive_counts(infolder, targettypelist, nontargetlist, tot_dict):
    genecounts = {}
    tallydict = {}

    poscountlines = []
    targetindices = []
    ntindices = []

    adjfactor = 0

    totalline = ['TOTAL', tot_dict['Target Tissues'],  tot_dict['ALL Non Target'], (tot_dict['Target Tissues'] + tot_dict['ALL Non Target'])]
    poscountlines.append(totalline)
    print(targettypelist)

    with open(f'{infolder}/PositiveCounts_CellType.csv') as poscountinfile:
        INposreader = csv.reader(poscountinfile)
        
        linecount = 0
        for line in INposreader:
            if linecount == 0:
                indice = 0

                for item in line:
                
                    if ':' in item:
                        letsmod = item.split(' : ')
                        upper_tiss_type = letsmod[0].capitalize()
                        
                        if "_" in upper_tiss_type:
                            upper_tiss_type = upper_tiss_type.replace('_', ' ')
                            
                        upper_cell_type = letsmod[1].capitalize()
                        actitemname = upper_tiss_type + ' : ' + upper_cell_type
                        
                    else:
                        actitemname = 'NULL'
                    
                    if actitemname in targettypelist:
                        targetindices.append(indice)
                    
                    elif indice > 0:
                        ntindices.append(indice)

                    indice += 1

            else:
                gene = line[0]

                genecounts[gene] = {}
                tallydict[gene] = [gene,]

                targetcount = 0
                ntcount = 0

                gene = line[0]

                adjfactor += 1

                for index in targetindices:
                    targetcount += int(line[index])

                genecounts[gene]['t'] = targetcount

                ubiquitytally = (targetcount/tot_dict['Target Tissues'])
                tallydict[gene].append(ubiquitytally)
                targetpercent = ubiquitytally * 100

                for index in ntindices:
                    ntcount += int(line[index])

                #Do the same for non target cells
                genecounts[gene]['nt'] = ntcount
                
                nontargettally = (ntcount/tot_dict['ALL Non Target'])
                tallydict[gene].append(1- (nontargettally))
                nontargetpercent = nontargettally * 100

                #Then we get a percent of ALL cells which express a given gene
                tot_positive = targetcount + ntcount
                genecounts[gene]['tot'] = tot_positive

                specificitytally = (targetcount/tot_positive)
                tallydict[gene].append(specificitytally)

                percentaretarget = (specificitytally)*100
                percentarenontarget = (ntcount/tot_positive)*100

                indivposcount = [gene, targetcount, ntcount, tot_positive, ' ', targetpercent, nontargetpercent,' ', percentaretarget, percentarenontarget]
                poscountlines.append(indivposcount)
        
            linecount += 1
    
    return tallydict, genecounts, poscountlines, adjfactor

#Makes a file containing the number of cells associated with each category (target, nontarget and total) which express each gene, 
#as well as Percerntages of each category which express the gene, and percentages of cells expressing the gene which belong to eahc category
def make_percentfile(outname, poscountlines):
    with open(f'{outname}/PercentPositiveCells_ByUserCategory.csv', 'w+') as poscountfile:
        OUTposcount = csv.writer(poscountfile)

        header = ['Gene', ' TARGET count', 'NONTARGET count', 'TOTAL count', ' ', '% Target Cells Expressing Gene', '% ALL Non Target Cells Expressing Gene', ' ', 'Percent cells expressing gene are Target Cells', 'Percent cells expressing gene are ANY NON-TARGET cell type']
        OUTposcount.writerow(header)

        for line in poscountlines:
            OUTposcount.writerow(line)

#Carries out a Chi-Squared test, with Bonferroni corrections, to check whether each gene is enriched  in target tissues relative to non target tissues. 
#This uses a STRICT pre-adjustment p-value threshold of p<0.01
def chisquare(genecounts, adjfactor, freq_t, freq_nt, tallydict):
    enrichmentdict = {}
    printables = []

    for gene in genecounts:
        targetcount = genecounts[gene]['t']
        ntcount = genecounts[gene]['nt']
        genetot = genecounts[gene]['tot']

        expcount_t = genetot * freq_t
        expcount_nt = genetot * freq_nt

        enrichmentdict[gene] = [[targetcount, ntcount], [expcount_t, expcount_nt]]
    
    pthresh = 0.01 / adjfactor

    for gene in enrichmentdict:
        actual, expected = enrichmentdict[gene]

        chsqresults = sp.chisquare(actual, f_exp = expected)
        stat = chsqresults[1]
        tally = 0
        wrongway = ''

        if stat < pthresh:
            if actual[0] < expected[0]:
                wrongway = 'YES'
            else:
                tally = 1
        
        outrow = [gene, tally, ' ', stat, wrongway]
        printables.append(outrow)

        tallydict[gene].append(tally)
    
    return tallydict, printables

#Produces a summary .csv containing the full statistical output of chi-squared tests
def print_stats(printables, adjfactor, freq_t, freq_nt, outname):
    with open(f'{outname}/EnrichmentInTargets_CHISQUAREDresults.csv', 'w+') as distribfile:
        pthresh = 0.01 / adjfactor

        OUTdistribution = csv.writer(distribfile)

        OUTdistribution.writerow([f'{adjfactor} chi-squared tests carried out'])
        OUTdistribution.writerow([f'Bonferroni adjusted p-value threshold of {pthresh} (Unadjusted is 0.01)'])
        OUTdistribution.writerow([f'Target Frequency = {freq_t}, Non-Target Frequency = {freq_nt}'])

        OUTdistribution.writerow([''])

        headrow = ['Gene Name', 'FCA_Significant Enrichment in Target Cells', ' ', 'p-value', 'Wrong direction?']
        OUTdistribution.writerow(headrow)

        for iline in printables:
            OUTdistribution.writerow(iline)

#desc
def handle_overlap_counts(infolder, genecounts, hittally):
    for gene in genecounts:

        allcount = genecounts[gene]['tot']

        overlaptally = 0
        numrefs = 0

        for file in glob.glob(f'{infolder}/OverlapCounts/OverlapCounts_{gene}.csv'):

            with open(file) as overlapsfile:
                INoverlaps = csv.reader(overlapsfile)
                next(INoverlaps)

                for line in INoverlaps:         
                    for x in line[1:]:
                        overlaptally += (float(x)/allcount)
                    
                    numrefs += 1

        finaloverlap = overlaptally/numrefs
        hittally[gene].append(finaloverlap)
    
    return hittally

#prints the final output containing FlyCellAtlas DIGITtally results
def print_output(outname, hittally, header):
    with open(f'{outname}/FlyCellAtlas_tally.csv', 'w+') as outputfile:
        OUTwriter = csv.writer(outputfile)

        OUTwriter.writerow(header)

        for entry in hittally:
            OUTwriter.writerow(hittally[entry])


def build_the_tally(fca_file_folder, celltypes, synonyms, outputatlas, outputtally, referencemode, associatedfiles):

    targets, allothers, totdict = get_celltypes(celltypes, associatedfiles)
    freq_t, freq_nt = get_freqs(totdict)

    tallydict, genecounts, poscountlines, adjfactor = handle_positive_counts(fca_file_folder, targets, allothers, totdict)
    make_percentfile(outputatlas, poscountlines)

    tallydict, chisquare_outs = chisquare(genecounts, adjfactor, freq_t, freq_nt, tallydict)
    print_stats(chisquare_outs, adjfactor, freq_t, freq_nt, outputatlas)

    header = ['Gene Name', 'Ubiquity Score', 'Proportion of Non-Target cells expressing gene', 'Proportion of cells expressing gene are Target cells', 'Significantly Enriched in Target Cells?']
    
    if referencemode == True:
        header.append('Overlap Score')

        tallydict = handle_overlap_counts(fca_file_folder, genecounts, tallydict)

    print_output(outputtally, tallydict, header)   
   
