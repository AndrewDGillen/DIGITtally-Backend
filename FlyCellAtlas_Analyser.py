#Code written by Dr. Andrew Gillen (Dow-Davies lab, University of Glasgow)
#Originally written October 2022, last updated 26/06/2023

import csv
import gc

import loompy as lp
import pandas as pd

#Part of the pipeline for DIGITtally analysis - see www.digittally.org
#Specifically, this program analyses the Single cell data in (https://flycellatlas.org/)

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
    
    #We gather the cell type associated with each cell
    for entry in cellcharacteristics:
        info = cellcharacteristics[entry]
        compiled = info[0] + ' : ' + info[1]
        if compiled not in unique_subsets:
            unique_subsets.append(compiled)
    
    #We then compile a dictionary counting the number associated with each cell type
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

	    #For each individual cell, we gather expression data seperately
            for cellid in cellcharacteristics:
                celltype = cellcharacteristics[cellid]
                compiledtype = celltype[0] + ' : ' + celltype[1]

                if compiledtype not in expression_by_subset:
                    expression_by_subset[compiledtype] = []
                    present_by_subset[compiledtype] = []

                goiread = goiframe.at[cellid,gene]

                expression_by_subset[compiledtype] += (goiread)

                #UMI >= 1 is used to classify a gene as "present" in FlyCellAtlas - see original paper (Hongjie et al., 2022)
                if goiread >= 1:
                
                    #If the gene is expressed in a given cell, we incrememnt the number of cells of that type by 1
                    present_by_subset[compiledtype] += 1
                    
                    #If reference genes are also in use, we check whether each reference gene is also expressed
                    if reflist != []:
                        for ref in expressiondict[cellid]:
                        
                            if ref not in overlap_by_subset:
                                overlap_by_subset[ref] = {}
                        
                            if compiledtype not in overlap_by_subset[ref]:
                                overlap_by_subset[ref][compiledtype] = 0
                            
                            #If both target and reference genes, we increment the count for that cell type by 1
                            if expressiondict[cellid][ref] >= 1:
                                overlap_by_subset[ref][compiledtype] += 1
                                
            #Now that overlap counts have been made, we check the % of each cell type which express each gene and each reference gene
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
            
            #We clear the dictionaries once useful data has been analysed to free up RAM
            del overlap_by_subset
            del goiframe
            gc.collect()
            
            goiframe = pd.DataFrame()
            overlap_by_subset = {}
            
            #We print the expression data from each gene to the output file
            countrow = [gene,]
            percentrow = [gene,]
            avgrow = [gene,]
            
            for subset in expression_by_subset:
            
                #If this is the first gene, we assemble the header listing cell types
                if topheaderprinted == 0:
                    toplevelheader.append(subset)
                
                #The total count by cell type is printed into the COUNT file
                totalcount = present_by_subset[subset]
                countrow.append(totalcount)
		
		#The % of cells of each type expressing a gene is printed to the PERCENTAGE file
                percent = (totalcount / subset_totals[subset]) * 100
                percentrow.append(percent)

		#The average expression of each gene in each cell type is printed to the AVERAGE file
                avgexpression = (expression_by_subset[subset]/subset_totals[subset])
                avgrow.append(avgexpression)
	    
	    #If this is the first gene, we print the header
            if topheaderprinted == 0:
                OUTcount.writerow(toplevelheader)
                OUTpercent.writerow(toplevelheader)
                OUTavg.writerow(toplevelheader)

                topheaderprinted += 1
            
            #We then print the gene's data
            OUTcount.writerow(countrow)
            OUTpercent.writerow(percentrow)
            OUTavg.writerow(avgrow)
	   
	    #Again, we free up RAM and reset dictionaries
            del present_by_subset
            del expression_by_subset
            gc.collect()
            
            present_by_subset = {}
            expression_by_subset = {}

def Analyse_FCA(loomfile, goi_checked, reference_gene_file, output_folder):

    #We gather lists of genes of interest and reference genes
    goilist, reflist = populate_lists(goi_checked, reference_gene_file)
 
    #The LOOM file is loaded
    df = lp.connect(loomfile, validate = False)
    
    #We gather cell/tissue/annotation types for each cell
    cellcharacteristics, subset_totals, unique_subsets = populate_singlecells(df)
    
    #Initialise output files
    make_total_out(subset_totals, output_folder)
    
    #We gather expression data for reference genes
    #IE - which cells express each gene
    expressiondict = find_ref_expression(reflist, cellcharacteristics, df)
    
    #Finally, we check the expression of each gene in each type of cell
    #Specifically, we check:	TOTAL NUMBER of cells expressing that gene; 	PERCENTAGE of that type that express that gene;		AVERAGE gene expression in that cell type
    find_goi_expression(goilist, expressiondict, cellcharacteristics, unique_subsets, subset_totals, df, output_folder, reflist)
    
    df.close()

