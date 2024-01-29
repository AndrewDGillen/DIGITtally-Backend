
import csv
import gc

import loompy as lp
import pandas as pd

#Creates lists of genes of interest and reference genes
def populate_lists(gois, refs):

    goilist = []
    reflist = []

    with open(gois) as goifile:
        for line in goifile:
            goilist.append(line[:-1])
    
    if refs != '':
        with open(refs) as reffile:
            for line in reffile:
                reflist.append(line[:-1])
    
    return goilist, reflist

#Creates necessary dictionaries from the .loom file to annotate incoming expression data
def populate_singlecells(loomdf):
    cellcharacteristics = {}
    subset_totals = {}

    unique_subsets = []

    cellarray = loomdf.ca["CellID", "tissue"]

    for item in cellarray:
        cellcharacteristics[item[0]] = [item[1],]

    cellarray = loomdf.ca["CellID", "annotation"]

    for item in cellarray:
        cellcharacteristics[item[0]].append(item[1])
    
    for entry in cellcharacteristics:
        info = cellcharacteristics[entry]
        compiled = info[0] + ' : ' + info[1]
        if compiled not in unique_subsets:
            unique_subsets.append(compiled)

    for subset in unique_subsets:
        subset_totals[subset] = 0

    for cellid in cellcharacteristics:
        celltype = cellcharacteristics[cellid]
        compiledtype = celltype[0] + ' : ' + celltype[1]

        subset_totals[compiledtype] += 1
    
    return cellcharacteristics, subset_totals, unique_subsets

#Very simple function to create a .csv containing the counts associated with each FlyCellAtlas cell type
def make_total_out(subset_totals, output):
    with open(f'{output}/Subtype_TotalCounts.csv', 'w+') as subtypetotalfile:
        OUTtotals = csv.writer(subtypetotalfile)
        totheader = ['Subset', 'Total Cells']
        OUTtotals.writerow(totheader)

        for subset in subset_totals:
            totoutrow = [subset, subset_totals[subset]]
            OUTtotals.writerow(totoutrow)

#Creates a dictionary containing the expression data for each reference gene in each cell - This is the limiting factor in reference gene list size due to the need to store values on EACH CELL in a dictionary
def find_ref_expression(reflist, cellcharacteristics, df):
    expressiondict = {}

    for reference in reflist:
        #This except will never be hit if the gene lists come from SCOPEscope, but it's maintained to allow the program to be used by itself, if that is desirable to anyone
        try:
            #A Pandas dataframe is created to hold the data from the loom file for each gene against cell ID
            refframe = pd.DataFrame(df[df.ra.Gene == reference,:], index=[reference], columns=df.ca["CellID"]).T
        except:
            print(f'\n{reference} not recognised within the LOOM file - have you run SCOPEscope?\n\nExiting until this issue is resolved.')
            exit()
        
        #Read data for EACH CELL is stored in expressiondict
        for cellid in cellcharacteristics:
            
            if cellid not in expressiondict:
                expressiondict[cellid] = {}

            read = refframe.at[cellid,reference]
            expressiondict[cellid][reference] = read

        #We clean up dataframes at each stage of the loop to try and decrease RAM requirement
        del refframe
        gc.collect()
        refframe = pd.DataFrame()

    return expressiondict

#The centrepiece function. Gathers expression data and defines each gene as present/absent in EACH CELL. This is then compared to reference gene presence.
def find_goi_expression(goilist, expressiondict, cellcharacteristics, unique_subsets, subset_totals, df, output, reflist):
    overlap_by_subset = {}

    #It seems pre-emptive to open these here but it allows us to create a single file for Average expression, Percent expression and Positive Counts without creating three huge dictionaries
    with open(f'{output}/AvgExpression_CellType.csv', 'w+') as avgbygenefile, open(f'{output}/PercentExpression_CellType.csv', 'w+') as percentbygenefile,open(f'{output}/PositiveCounts_CellType.csv', 'w+') as countbygenefile:
        OUTavg = csv.writer(avgbygenefile)
        OUTpercent = csv.writer(percentbygenefile)
        OUTcount = csv.writer(countbygenefile)

        toplevelheader = ['Gene',]
        topheaderprinted = 0
        
        #Every gene of interest is processed seperately
        for gene in goilist:
            
            expression_by_subset = {}
            present_by_subset = {}

            header = ['Reference Gene vs Cell Type']

            for subset in unique_subsets:
                expression_by_subset[subset] = 0
                present_by_subset[subset] = 0

            try:
                goiframe = pd.DataFrame(df[df.ra.Gene == gene,:], index=[gene], columns=df.ca["CellID"]).T
            except:
                print(f'\n{reference} not recognised within the LOOM file - have you run SCOPEscope?\n\nExiting until this issue is resolved.')
                exit()


            for cellid in cellcharacteristics:
                celltype = cellcharacteristics[cellid]
                compiledtype = celltype[0] + ' : ' + celltype[1]

                if compiledtype not in expression_by_subset:
                    expression_by_subset[compiledtype] = []
                    present_by_subset[compiledtype] = []

                goiread = goiframe.at[cellid,gene]

                expression_by_subset[compiledtype] += (goiread)

                #UMI >= 1 is used to classify a gene as "present" in FlyCellAtlas - see original paper ()
                if goiread >= 1:
                    present_by_subset[compiledtype] += 1
                    if reflist != []:
                        for ref in expressiondict[cellid]:
                        
                            if ref not in overlap_by_subset:
                                overlap_by_subset[ref] = {}
                        
                            if compiledtype not in overlap_by_subset[ref]:
                                overlap_by_subset[ref][compiledtype] = 0
                        
                            if expressiondict[cellid][ref] >= 1:
                                overlap_by_subset[ref][compiledtype] += 1
            if reflist != []:
                with open(f'{output}/OverlapPercents/PercentOverlap_{gene}.csv', 'w+') as overlapbygenefile, open(f'{output}/OverlapCounts/OverlapCounts_{gene}.csv','w+') as overlapcountfile:
                    OUTpercentoverlapgene = csv.writer(overlapbygenefile)
                    OUTcountoverlapgene = csv.writer(overlapcountfile)

                    subheadersprinted = 0

                    for reference in overlap_by_subset:
                        countoutline = [reference]
                        percentoutline = [reference]
                
                        for celltype in overlap_by_subset[reference]:
                        
                            overlapping = overlap_by_subset[reference][celltype]
                            percentol = (overlapping/subset_totals[celltype]) * 100

                            if subheadersprinted == 0:
                                header.append(celltype)
                        
                            countoutline.append(overlapping)
                            percentoutline.append(percentol)
                    
                        if subheadersprinted == 0:
                            OUTpercentoverlapgene.writerow(header)
                            OUTcountoverlapgene.writerow(header)

                            subheadersprinted += 1
                    
                        OUTpercentoverlapgene.writerow(percentoutline)
                        OUTcountoverlapgene.writerow(countoutline)

            del overlap_by_subset
            del goiframe
            gc.collect()
            
            goiframe = pd.DataFrame()
            overlap_by_subset = {}

            countrow = [gene,]
            percentrow = [gene,]
            avgrow = [gene,]

            for subset in expression_by_subset:
                if topheaderprinted == 0:
                    toplevelheader.append(subset)

                totalcount = present_by_subset[subset]
                countrow.append(totalcount)

                percent = (totalcount / subset_totals[subset]) * 100
                percentrow.append(percent)

                avgexpression = (expression_by_subset[subset]/subset_totals[subset])
                avgrow.append(avgexpression)

            if topheaderprinted == 0:
                OUTcount.writerow(toplevelheader)
                OUTpercent.writerow(toplevelheader)
                OUTavg.writerow(toplevelheader)

                topheaderprinted += 1
            
            OUTcount.writerow(countrow)
            OUTpercent.writerow(percentrow)
            OUTavg.writerow(avgrow)

            del present_by_subset
            del expression_by_subset
            gc.collect()
            present_by_subset = {}
            expression_by_subset = {}

def Analyse_FCA(loomfile, goi_checked, reference_gene_file, output_folder):

    goilist, reflist = populate_lists(goi_checked, reference_gene_file)

    df = lp.connect(loomfile, validate = False)

    cellcharacteristics, subset_totals, unique_subsets = populate_singlecells(df)
    make_total_out(subset_totals, output_folder)

    expressiondict = find_ref_expression(reflist, cellcharacteristics, df)

    find_goi_expression(goilist, expressiondict, cellcharacteristics, unique_subsets, subset_totals, df, output_folder, reflist)
    
    df.close()

