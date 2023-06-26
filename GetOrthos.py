#Code written by Dr. Andrew Gillen (Dow-Davies lab, University of Glasgow)
#Originally written October 2022, last updated 26/06/2023

import argparse
import time
import csv

from collections import defaultdict
from tqdm import tqdm

from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.common.exceptions import TimeoutException

#Part of the pipeline for DIGITtally analysis - see www.digittally.org
#Not actually called during website running, but used to generate lists of orthologs for each Drosophila gene

#Interactions with websites are handled via the Web_Interface class. Functions are definted below.
class Web_Interface():
    def __init__(self, species, genelist, aliasdict, keggdict):
        self.genelist = genelist
        self.species = species
        self.aliasdict = aliasdict
        self.keggdict = keggdict

    #Creates a selenium "driver" object for flybase (standardised web page)
    def initialise_flybase_driver(self):
        options = Options()
        options.add_argument('-headless')
        fb_driver = webdriver.Firefox()
        fb_driver.get('https://flybase.org/')

        return fb_driver
    
    #Creates a selenium "driver" object for orthodb or kegg (individual links by gene)
    def initialise_other_driver(self, link):
        options = Options()
        options.add_argument('-headless')
        odb_driver = webdriver.Firefox()
        odb_driver.get(link)

        return odb_driver
    
    #A funtion to find the "orthologs" table on flybase. This is only run the first time (for each set of 10 genes processed).
    def find_first_table(self, driver):
        driver.find_element_by_css_selector(".col-xs-12[data-toggle-target='orthologs_sub").click()
        time.sleep(1)

        driver.find_element_by_css_selector(".col-xs-12[data-toggle-target='orthodb_orthologs_sub']").click()
        time.sleep(5)

    #Selenium occasionally runs into problems identifying the ortholog table on successive page loads.
    #To circumvent this, on sequential loads of the flybase page, we open and close the table tab
    def refresh_table(self, driver):
        expand_others = driver.find_element_by_css_selector(".col-xs-12[data-toggle-target='orthodb_orthologs_sub']")
        expand_others.click()
        expand_others.click()
        time.sleep(5)

    #Searches for a given gene on Flybase, then finds and opens the "other Species" orthology tab
    def find_ortho_table(self, driver, gene):

        #finds the search box on flybase and pastes our gene name in
        driver.find_element_by_id("GeneSearch").send_keys(gene)
        
        #clicks on the search button
        driver.find_element_by_id('GeneSearch_submit').click()
        
        time.sleep(3)

        #runs the appropriate function to locate and open/refresh the orthology tab
        if self.flybase_open == 0:
            self.find_first_table(driver)
        
        else:
            try:
                self.refresh_table(driver)

            except:
                #I have never triggered this except clause  - my guess is it shouldn't happen unless there are some weird genes deep in FlyBase without an Ortholog table 
                # or the FlyBase CSS is altered.
                print('ERROR: NOT FINDING TABLE')
                exit()

        try:
            WebDriverWait(driver, 5).until(EC.presence_of_element_located((By.CSS_SELECTOR, "table.mantine-Table-root")))
            orthotable  = driver.find_elements_by_css_selector("table.mantine-Table-root")

        #this clause runs when the ortholog table cannot be found in 5 seconds
        except TimeoutException:

            try:
                #notable isn't actually used - this is just a check to see whether the table is actually absent or the page has just failed to load
                notable = driver.find_element_by_css_selector(".mantine-Text-root")
                return "NO ORTHOLOGS", driver

            except:
                #In the occasion a page hasn't loaded, an error is returned so the webpage will be reloaded.
                raise TimeoutError
        
        return orthotable, driver

    #This function gathers useful data from the OrthoDB page for individual genes, if necessary
    def grab_useful(self, odb_driver, is_bombyx):
        
        #This function waits until the table of OrthoDB information can be found. 
        #If this times out, this function will exit with a TimeoutException, which is expected behaviour
        WebDriverWait(odb_driver, 5).until(EC.presence_of_element_located((By.CSS_SELECTOR, "div.s-group-ortho-annotations:nth-child(4)")))
        found_ortho_info = odb_driver.find_elements_by_css_selector('div.s-group-ortho-annotations:nth-child(4)')

        #each item of interest is defaulted to "null"
        entrezno = 'null'
        ensdata = 'null'
        keggno = 'null'

        #The orthoDB table is parsed, then each line is handled seperately
        for item in found_ortho_info:
            altdata = item.text
            ads = altdata.split('\n')
            
            #We check every line within the OrthoDB table to see if it's something we can make use of
            for subitem in ads:
                
                #The entrez number is gathered if present, in case it is needed for KEGG identification
                if 'Entrez:' in subitem:
                    entrezno = subitem.split(' ')[1]
                    entrezno = entrezno[3:]

                    try:
                        entrezno = entrezno.split(' ')[0]
                    except:
                        pass
                
                #For A. gam and A. aeg, the Ensembl ID is the sole thing we need - and indeed the only one we can use. This can also be used to find SilkDB identifiers
                if 'Ensembl' in subitem:
                    ensdata =  subitem.split(': ')[1]

                    try:
                        ensdata = ensdata.split(' ')[0]
                    except:
                        pass
                
                #OrthoDB can be fairly spotty with B. mori Ensembl IDs, and so we can sometimes use KEGG identifiers as an alternative route to identify orthologs
                if 'KEGGpathway' in subitem and is_bombyx == True:
                    keggpath = subitem.split(' ')[1]
                    keggpath = keggpath.split(';')[0]

                    #we use the first KEGG pathway associated with an ortholog as a means to derive its kegg identifier
                    buttons = odb_driver.find_elements_by_xpath(f"//*[contains(text(), '{keggpath}')]")
                    
                    #The link found on the OrthoDB page is used to navigate to KEGG and gather the kegg identifier
                    for button in buttons:
                        kegglink = button.get_attribute("href")

                        complete = 0

                        while complete == 0:
                            keggno, complete = self.get_kegg_info(kegglink, entrezno)
            
            return ensdata, keggno
    
    #This function gathers useful data from the KEGG page for individual genes, if necessary - THIS IS ONLY RUN FOR BOMBYX MORI ORTHOLOGY
    def get_kegg_info(self, link, entrezno):
        kegg_driver = self.initialise_other_driver(link)
        
        keggno = ''

        #Again, we wait to load the KEGG page fully, and if loading times out an exit condition of 0 is returned, prompting a retry
        try:
            WebDriverWait(kegg_driver, 10).until(EC.presence_of_element_located((By.CSS_SELECTOR, "body > div:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > table:nth-child(1) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(3)")))
            
            #If the supplied entrez number is found on the KEGG page, we locate all lines in which it is present (should only be one)
            if entrezno in kegg_driver.page_source:
                find_in_kegg = kegg_driver.find_elements_by_xpath(f"//*[contains(text(), '{entrezno}')]")

                #This just handles the specific layout of elements on KEGG - once the entrez id has been found, we find the associated KEGG identifier
                for item in find_in_kegg:
                    parentelement = item.find_element_by_xpath("..//..//..")
                    fullrow = parentelement.text
                    keggno = fullrow.split('KO:')[1]

                    if ' ' in keggno:
                        keggno = keggno.split(' ')[0]
                        
                    keggno = keggno[:-1]

            #If the page loads correctly, but the entrez id can't be found on the page, we just return an empty string, along with the "complete" exit condition.
            
            complete = 1
        
        except TimeoutException:
            complete = 0

        kegg_driver.close()

        return keggno, complete

    #Incredibly simple check for whether the species is Bombyx mori as, due to the structure of SilkDB, these must be processed differently
    def is_bombyx(self, species):
        if species == 'Bombyx mori':
            return True
        else:
            return False

    #Co-ordinates interface with other webpages (OrthoDB, KEGG) to follow up on genes whose orthologs cannot be found solely from FlyBase
    def followup(self, link, species):
        odb_driver = self.initialise_other_driver(link)

        #if grab_useful fails, an exit condition of 0 is returned, and the followup command will be repeated
        try:
            ensdata, keggno = self.grab_useful(odb_driver, self.is_bombyx(species))

        except TimeoutException:
            odb_driver.close()
            return [], 0

        odb_driver.close()

        #if neither enzembl id nor kegg id can be found for a given gene, "null" is returned, signifying no orthologs can be identified
        if ensdata == 'null' and keggno == 'null':
            return ensdata, 1

        #if either data source has been utilised to get a usable identifier, we process these into our desired format
        else:
            return self.process_found_data(species, ensdata, keggno)

    #Processes the found ortholog identifiers into a single "orthoname" list 
    def process_found_data(self, species, ensdata, keggno):
        orthoname = []

        #Sometimes, multiple ensembl ids are listed for a gene. in this case, these are split up and processed seperately
        if ';'  in ensdata:
            ensdata = ensdata.split(';')
        else:
            ensdata = [ensdata]

        #Each individual identifier in ensdata is processed seperately
        for entry in ensdata:

            #if the species is B. mori, we need to try to match up available info with a SilkDB identifier. 
            if self.is_bombyx(species):

                #First, we try to match up the ensembl ID
                try:
                    orthoname = orthoname + self.aliasdict[entry]

                #If this fails, we try the KEGG identifier
                except:
                    try:
                        orthoname = orthoname + self.keggdict[keggno]
                    
                    #If this fails, we're out of ways to match up with the SilkDB identifier, and so an empty orthoname is returned
                    except:
                        pass
            
            #For the other species, in which ensembl ID is directly usable, these are added straight to the Orthoname list
            else:
                orthoname.append(entry)

        return orthoname, 1

    #Searches for orthology information from each species of interest within the identified orthology table
    def query_table(self, species, orthotable):
        row_count = 1
        matchingorthos = []

        #This try/except is probably unneccessary, but accounts for cases in which the orthology table doesn't contain rows - 
        #IE contains no orthology information
        try:
            #For each row in the table, we check whether the species is one we're interested in
            for indivrow in orthotable[0].find_elements_by_css_selector("tr"):
                position = 1

                for cell in indivrow.find_elements_by_css_selector('td'):
                    #Species name is held in column 1
                    if position == 1:
                        speciesname = cell.text
                    
                    #ortholog ID is held in column 3
                    elif position == 3:
                        ortholog = cell.text
                        link = cell.find_element_by_xpath(".//a").get_attribute("href")
                    
                    position += 1
                
                #For every line past the header, we check if the species name matches our species of interest
                if row_count > 1:
                    if speciesname == species:
                        
                        #some ortholog IDs are lists - these are processed seperately
                        if ';' in ortholog:
                            ortho_list = ortholog.split(';')
                        else:
                            ortho_list = [ortholog]
                        
                        for indiv_ortho in ortho_list:
                            
                            #In A. gam/A. aeg cases where the ortholog ID does not contain the species identifier (ie is not in ensembl format), it is necessary to harvest this ensembl ID from OrthoDB
                            #Due to the nature of SilkDB and the way FlyBase displays B. mori ortholog IDs, it is ALWAYS necessary to consult OrthoDB for B.mori orthologs
                            if self.short_id not in indiv_ortho or species == 'Bombyx mori':
                                done = 0

                                #OrthoDB followup is attempted until successful
                                while done == 0:
                                    indiv_ortho, done = self.followup(link, species)

                                #If followup executes correctly, it should produce a list.
                                if type(indiv_ortho) == list:
                                    for indivmatch in indiv_ortho:

                                        #as long as there's any identifier provided for the ortholog, it is added to matchingorthos
                                        if indivmatch != '':
                                            matchingorthos.append(indivmatch)

                                #I don't actually know if this runs tbh
                                else:

                                    if indiv_ortho not in matchingorthos and indiv_ortho not in  ['null', 'No Orthology information available'] :
                                        print(f'\n\n **** THIS BLOCK IS TRIGGERED ****\n\n')
                                        print(matchingorthos)
                                        matchingorthos.append(indiv_ortho)
                            
                            #If the ensembl id is found, and this is enough to use, we just add this to matchingorthos
                            else:
                                matchingorthos.append(indiv_ortho)

                row_count += 1

        #If orthology information CANNOT BE FOUND - ie no table on flybase, we report this
        except:
            matchingorthos.append('No Orthology information available')
        
        #On the other hand, if the table EXISTS, but NO ORTHOLOGS ARE FOUND IN A GIVEN SPECIES, the returned message differs slightly
        if matchingorthos == []:
            matchingorthos.append('No orthologs found in this species')
        
        return matchingorthos

    #This function drives the web search, co-ordination ortholog finding
    def get_orthos(self):

        #This dictionary is manually created based on the short species identifiers used in entrez IDs
        short_ids = {
            'Aedes aegypti' : 'AAEL',
            'Anopheles gambiae' : 'AGAP',
            'Bombyx mori' : 'BGIBMGA'
        }

        # A dictionary of dictionaries is created to hold orthology data in the format: {Species:{gene:ortholog}}
        match_dict = defaultdict(dict)

        self.flybase_open = 0
        
        #A webdriver is created for FlyBase
        fb_driver = self.initialise_flybase_driver()

        #Data is gathered for all species for one gene at a time, to minimise reloading of FlyBase
        for gene in self.genelist:
            
            #A while loop is used to repeat the orthology table finding process until either it is found or absence is confirmed by presence of the "No orthology data available" message.
            found = 0
            while found == 0:
                try:
                    ortho_table, fb_driver = self.find_ortho_table(fb_driver, gene)
                    found = 1

                except:
                    pass

            self.flybase_open += 1

            #Each species is handled seperately
            for species in self.species:
                self.short_id = short_ids[species]

                #As long as the orthology table is found, it will be queried for orthologs in each species of interest
                if ortho_table != "NO ORTHOLOGS":
                    match_dict[species][gene] = set(self.query_table(species, ortho_table))

                #Otherwise, the gene is listed as possessing no orthology information
                else:
                    match_dict[species][gene] = {'No Orthology information available'}

                ###THIS IS FOR TESTING###
                print(' ')
                print(gene,  species, match_dict[species][gene])
                print(' ')

        fb_driver.close()
        return match_dict

#Simple function for parsing genes of interest
def populate_lists(annos):

    annolist = []

    with open(annos) as annofile:
        for line in annofile:
            annolist.append(line.strip('\n'))
    
    return annolist

#SilkDB, for whatever reason, uses its own identifier for genes rather than ensembl ids. 
#To work around this, the ensembl IDs and KEGG identifiers associated with each SilkDB identifier are held in dictionaries
def pop_silkdb_info():
    aliasdict = {}
    keggdict = {}

    with open('/media/sf_SHARED_FOLDER/silkdb_annotation.txt') as silkdbfile:
        
        for line in silkdbfile:
            expandedinfo = line.split('\t')

            silkdbname = expandedinfo[0]

            alias = expandedinfo[3].split(',')
            alias = alias[0]

            aliasdict[alias] = silkdbname

            kegg = expandedinfo[7]

            if kegg not in keggdict:
                keggdict[kegg] = [silkdbname,]

            else:
                keggdict[kegg].append(silkdbname)
    
    return aliasdict, keggdict

#Gathers user arguments
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', type = str, default = 'TGF/Genelists/DetectedGOIs_FbID.txt', help = 'gene list (.txt) in FlyBase ID format (Gene Symbols may also work though this is less consistent)')
    parser.add_argument('-s', '--species', nargs = '+', default = ['Aedes aegypti', 'Anopheles gambiae', 'Bombyx mori'], choices = ['Aedes aegypti', 'Anopheles gambiae', 'Bombyx mori'], help = 'Species to get orthologs for')
    parser.add_argument('-m', '--mothfile', type = str, default = 'TGF/AssociatedFiles/SilkDB_Annotations.txt')
    parser.add_argument('-o', '--output', type = str, default = 'TGF/Orthology', help = 'Give desired Orthology output folder name.')
    return parser.parse_args()

#creates ortholog output files for each species and prints the header row. 
def create_files(species, outfolder):

    for individual_species in species:

        with open(f'{outfolder}/{individual_species}_orthologs.csv', 'w+') as orthologfile:
            OUTorthos = csv.writer(orthologfile)

            header = ['Gene Name', 'Orthologs']
            OUTorthos.writerow(header)

#Amends the ortholog files with discovered information
def update_file(species, outfolder, matchingorthos):
    with open (f'{outfolder}/{species}_orthologs.csv', 'a+') as orthologfile:
            OUTorthos = csv.writer(orthologfile)

            for indivgene in matchingorthos:

                if len(matchingorthos[indivgene]) > 1:
                    orthos_to_print = ", ".join(matchingorthos[indivgene])
                else:
                    orthos_to_print = ''.join(matchingorthos[indivgene])

                OUTorthos.writerow([indivgene, orthos_to_print])

#carries out searches for a given list of genes and appends discovered orthology information to the orthology out files
def execute_search(species, genelist, aliasdict, keggdict, outfolder):

    orthos = Web_Interface(species, genelist, aliasdict, keggdict)
    matchingorthos = orthos.get_orthos()

    for individual_species in species:
        update_file(individual_species, outfolder, dict(matchingorthos[individual_species]))


def main():

    args = get_args()
    start_time = time.time()  
    
    #Creates output files
    create_files(args.species, args.output)
    
    #silkdb info is only populated if needed
    if 'Bombyx mori' in args.species:
        aliasdict, keggdict = pop_silkdb_info()
    else:
        aliasdict, keggdict = {}, {}

    #The genelist is always populated
    genelist = populate_lists(args.input)

    tempgenelist = []
    totalprocessed = 0

    #Genes are processed TEN AT A TIME just to ensure selenium crashing, which is frustratingly common, doesn't wipe ALL progress.
    for gene in tqdm(genelist):
        tempgenelist.append(gene)
        totalprocessed += 1

        if len(tempgenelist) == 10 or totalprocessed == len(genelist):        
            execute_search(args.species, tempgenelist, aliasdict, keggdict,args.output)
            tempgenelist = []

    print('--- %s seconds ---' % (time.time() - start_time))
  
main()
