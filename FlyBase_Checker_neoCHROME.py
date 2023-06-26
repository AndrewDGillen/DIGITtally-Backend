#Code written by Dr. Andrew Gillen (Dow-Davies lab, University of Glasgow)
#Originally written October 2022, last updated 26/06/2023

import csv
import time
import os
import shutil

from selenium import webdriver
from selenium.webdriver.support.ui import Select
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.proxy import *
from fake_useragent import UserAgent

#Part of the pipeline for DIGITtally analysis - see www.digittally.org
#Specifically, this program analyses the literature data in FlyBase (http://flybase.org/)

#Gathers CONFIG info
#Config file location is not listed here for security purposes
import json
with open("CONFIG FILE") as config_file:
    config = json.load(config_file)

#The interface with FlyBase is handled in its own class, utilising selenium and chromedriver for Chrome
#Unfortunately, this information is NOT QUERYABLE OFFLINE - thus there is an absolute requirement for internet connectivity via Selenium/Chrome
class Web_Interface():
    def __init__(self, tissue, outpath, writable_annotation):
        self.tissue = tissue
        self.filename = writable_annotation
        #ensures correct structure of flybase output
        self.FBdownload = f"{outpath}/Downloads"
        self.outpath = f"{outpath}/Downloaded_lists"
    
    #Uploads the genelist to FlyBase, selects synonyms, and downloads the output
    def find_annotated_genes(self):

        ua = UserAgent()
        userAgent = ua.chrome

        errors = [0, 0]
        
        chrome_options = Options()
        chromedriver = f'{config["FILE_LOCATION"]}/DIGITtallyScripts/chromedriver'

        chrome_options.add_argument('--no-sandbox')
        chrome_options.add_argument("--disable-extensions") 
        chrome_options.add_argument("--disable-gpu") 
        chrome_options.add_argument('--disable-dev-shm-usage')
        chrome_options.add_argument('--ignore-certificate-errors')


        chrome_options.add_argument("--disable-setuid-sandbox")
        chrome_options.add_argument('--allow-insecure-localhost')
        prefs = {
        'download.default_directory' : self.FBdownload,
        'credentials_enable_service': False,
        "profile.default_content_setting_values.notifications" : 2,
        'profile': {'password_manager_enabled': False}
        }
        chrome_options.add_argument(f'user-agent={userAgent}') 
        chrome_options.add_experimental_option('prefs', prefs)

        chrome_options.add_argument("window-size=1920,1080")
        
        from selenium.webdriver.common.desired_capabilities import DesiredCapabilities
        polipo_proxy = "http://wwwcache.gla.ac.uk:8080"
        proxy = Proxy({
            'proxyType': ProxyType.MANUAL,
            'httpProxy': polipo_proxy,
            'httpsProxy':polipo_proxy,
            'ftpProxy' : polipo_proxy,
            'sslProxy' : polipo_proxy,
            'noProxy'  : ''
        })
        desired_capabilities = dict(DesiredCapabilities.CHROME)
        proxy.add_to_capabilities(desired_capabilities)
        desired_capabilities['acceptInsecureCerts'] = True         
        desired_capabilities['acceptSslCerts'] = True

        chrome_options.add_argument('--disable-blink-features=AutomationControlled')
        chrome_options.add_experimental_option("excludeSwitches", ["enable-automation"])
        chrome_options.add_experimental_option('useAutomationExtension', False)
        chrome_options.add_argument("start-maximized")
        chrome_options.add_argument('--headless')
        driver = webdriver.Chrome(executable_path=chromedriver, options=chrome_options, desired_capabilities=desired_capabilities)
   

        #driver.execute_script("Object.defineProperty(navigator, 'webdriver', {get: () => undefined})")
        #driver.execute_cdp_cmd('Network.setUserAgentOverride', {"userAgent": 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/83.0.4103.53 Safari/537.36'})
        print(driver.execute_script("return navigator.userAgent;"))
              
        driver.get('https://flybase.org/cgi-bin/qb.pl')
        time.sleep(30)
        html_source = driver.page_source
        jscript_source = driver.execute_script("return document;")
        choose_qb = driver.find_element("css selector", '#segment01 > b:nth-child(2)')

        choose_qb.click()

        time.sleep(10)

        #Tells flybase we want to use the Controlled Vocabulary option
        use_cvs = Select(driver.find_element("css selector",'#editorFull > table:nth-child(5) > tbody:nth-child(1) > tr:nth-child(1) > td:nth-child(1) > select:nth-child(2)'))
        use_cvs.select_by_value('go')
        
        time.sleep(10)

        #Sends the current tissue annotation to the search box
        desired_anno = driver.find_element(By.CSS_SELECTOR, '#editorFull > table:nth-child(10) > tbody:nth-child(1) > tr:nth-child(3) > td:nth-child(1) > input:nth-child(2)')
        desired_anno.send_keys(self.tissue)

        startsearch = driver.find_element(By.CSS_SELECTOR, '#searchhier')
        startsearch.click()

        time.sleep(10)
	
	#Finds the correct location for the CV term selection
        correct_row = driver.find_element("css selector",'#editorFull > table:nth-child(10) > tbody:nth-child(1) > tr:nth-child(3) > td:nth-child(1)')

        test = driver.find_element(By.CSS_SELECTOR, '#editorFull > table:nth-child(10) > tbody:nth-child(1) > tr:nth-child(3) > td:nth-child(1) > input:nth-child(5)')

        input_set = driver.find_element(By.CSS_SELECTOR,'#editorFull > table:nth-child(10) > tbody:nth-child(1) > tr:nth-child(3) > td:nth-child(1)')

        for element in input_set.find_elements(By.CSS_SELECTOR, '*'):


            if element.get_attribute('value') != None:
                if element.get_attribute('value') == "Use this term ID":
                    current_link = 'input.lightgreen:nth-child(4)'

                if ';' in  element.get_attribute('value'):
                    split_element = element.get_attribute('value').split(' ; ')

                    if split_element == desired_anno:
                        target_to_click = current_link
                        break
        
        print(split_element, current_link)

        correct_search_func = driver.find_element("css selector", current_link)
        correct_search_func.click()

        time.sleep(10)

        #Runs the search
        run = 0
        tries = 0
        while run == 0:
            try:
                beginsearch = driver.find_element("xpath",'//*[@id="run"]')
                run = 1
            except:
                if tries > 5:
                    exit()
                print("can't find")
                time.sleep(5)
                tries += 1
             

        beginsearch.click()

        time.sleep(10)

        #As this program is set up to use the FlyBase Controlled Vocabulary, there should NEVER be an issue where a search term is not recognised
        #AS LONG AS THE ONTOLOGY FILE IS CORRECT - PLEASE USE TGF_EX_CreateFlyBaseOntologies or the prebuilt FlyBase_Ontologies.csv
        basepage = driver.current_window_handle
        oldtabs = driver.window_handles
        
        #The search finds two sets of pages that we're interested in - Genes with the annotation on the page anywhere ("Any known link to...")
        #And alleles associated with a phenotype in the annotated type ("Any allele causing a phenotype associated with...")

        try:
            do_genes = driver.find_element("css selector",'#alt-fbgn')
            print(driver.current_url)
            print(driver.execute_script("return navigator.webdriver;"))
            do_genes.click()
            print(driver.current_url)
            time.sleep(60)
            #These links open in new tabs, so we keep the old page open and switch to each in turn.
            newtabs = driver.window_handles
            
            print(oldtabs)
            print(newtabs)
            
            for tab in newtabs:
                if tab in oldtabs:
                    pass
                else:
                    print(tab)
                    print(driver.execute_script("return navigator.webdriver;"))
                    new_tab=tab

            time.sleep(20)        
            driver.switch_to.window(new_tab)
            time.sleep(30)


            #The list of genes which contain the annotation on their page is downloaded using the dropdown message
            download_genes = driver.find_element("css selector",'#hitlist-export')
            
            download_genes.click()
            
            time.sleep(2)
            goahead = driver.find_element("css selector",'ul.dropdown-menu:nth-child(3) > li:nth-child(8) > a:nth-child(1)')
            
            goahead.click()
            
            time.sleep(10)

            #Moves the downloaded hitlist to the output folder, then closes the window
            complete = False

            while complete == False:
                try:
                    shutil.move(f'{self.FBdownload}/FlyBase_IDs.txt', f'{self.outpath}/{self.filename}_AnyAnnotations.txt')
                    complete = True
                    print(self.tissue, 'genes')
                except:
                    time.sleep(10)

            driver.close()
            driver.switch_to.window(basepage)
        
        #THis handles cases where no genes are possess a given annotation. Why are these annotations in the FlyBase ontology then? No one knows.
        except:
            try:
                
                test_point = driver.find_element(By.CSS_SELECTOR, "#wrapper_two > table:nth-child(4)")

                print(f'No genes associated with tissue {self.tissue}')

                errors[0] = 1
            
            except Exception as e:

                print(e)
                raise LookupError
            
        print(driver.current_url)
	
        #The same process is repeated for the alleles associated with the annotation
        try:
            do_genes = driver.find_element("css selector",'#alt-fbal')
            do_genes.click()
            
            oldtabs = newtabs
            
            newtabs = driver.window_handles
            print(newtabs)
            for tab in newtabs:
                if tab in oldtabs:
                    pass
                else:
                    print(driver.execute_script("return navigator.webdriver;"))
                    new_tab=tab
                    break

            driver.switch_to.window(new_tab)
            print(driver.current_url)
            
            time.sleep(5)

            download_genes = driver.find_element("css selector",'#hitlist-export')
            download_genes.click()
            
            time.sleep(2)
            goahead = driver.find_element("css selector",'ul.dropdown-menu:nth-child(3) > li:nth-child(6) > a:nth-child(1)')
            goahead.click()

            time.sleep(10)

            complete = False
            while complete == False:
                try:
                    shutil.move(f'{self.FBdownload}/FlyBase_IDs.txt', f'{self.outpath}/{self.filename}_AllelesCausingPhenotypes.txt')
                    complete = True
                    print(self.tissue, 'alleles')
                except:
                    time.sleep(10)
            driver.close()
            
            driver.switch_to.window(basepage)
            print(driver.current_url)

        #THis handles cases where no alleles are associated with a phenotype in a given annotation
        except:
            try:
                
                test_point = driver.find_element(By.CSS_SELECTOR, "#wrapper_two > table:nth-child(4)")

                print(f'No alleles associated with tissue {self.tissue}')

                errors[1] = 1

            except Exception as e:

                print(e)
                raise LookupError
        
        driver.close()

        return errors

#Creates a list of the desired annotations, from the annotation text file.
def populate_lists(annos):

    annolist = []

    with open(annos) as annofile:
        for line in annofile:
            annolist.append(line.strip('\n'))
    
    return annolist

#Uses the supplied FlyBase Alleles ⇔ Genes file to associate each allele found 
def convert_alleles(alleles, allelefile):
    converted = []

    with open(allelefile) as af:
        INallele = csv.reader(af, delimiter = '\t')
        
        next(INallele)
        next(INallele)

        for line in INallele:
            fbal = line[0]
            fbid = line[2]

            if fbal in alleles:
                converted.append(fbid)
    
    return converted

#Parses through downloaded lists to find whch genes from the gene of interest list:
	#Are associated with a given phenotype AT ALL
	#have an allele which is specifically causative of a specific phenotype
def check_lists(annotation, outpath, fulloutdict, tallyany, tallyphenotype, genelist, allelefile, error_cases, workinganno):

    fulloutdict['Header'].append(f'Any known link to {annotation}?')
    fulloutdict['Header'].append(f'Any allele causing a phenotype in {annotation}?')
    
    time.sleep(5)
    #Parses "ANY ASSOCIATION" data
    if error_cases[0] == 0:
        print(f'{outpath}/Downloaded_lists/{workinganno}_AnyAnnotations.txt')
        with open(f'{outpath}/Downloaded_lists/{workinganno}_AnyAnnotations.txt') as file1:
            anygenes = file1.readlines()
    else:
        print('cannot open file')
        anygenes = []
    
    #Parses "PHENOTYPE" data
    if error_cases[1] == 0:
        print(f'{outpath}/Downloaded_lists/{workinganno}_AllelesCausingPhenotypes.txt')
        with open(f'{outpath}/Downloaded_lists/{workinganno}_AllelesCausingPhenotypes.txt') as file2:
            phenogenes = file2.readlines()
    else:
        print('cannot open file')
        phenogenes = []
    
    #We check each gene in the gene of interest list against the list of genes associated with a term
    #If a gene is associated with a term, it scores a point. If not, it doesn't
    found = []
    alleles = []
    for line in anygenes:
        geneid = line.strip('\n')
        
        if geneid in genelist:
            tallyany[geneid] += 1
            fulloutdict[geneid].append(1)
            found.append(geneid)

    for gene in genelist:
        if gene not in found:
            fulloutdict[gene].append(0)
    print('checking of "Any Genes" completed')

    
    #Next, we create a list of alleles causing alterations in a given term
    #These alleles are connverted to a list of genes, then this list is checked against genes of interest
    #If a gene is associated with a term, it scores a point. If not, it doesn't
    found2 = []
    for line in phenogenes:
        alleleid = line.strip('\n')
        alleles.append(alleleid)
    
    convertedgenes = convert_alleles(alleles, allelefile)

    for geneid in convertedgenes:
        if geneid in genelist and geneid not in found2:
            tallyphenotype[geneid] += 1
            fulloutdict[geneid].append(1)
            
            if fulloutdict[geneid][-2] == 0:
                tallyany[geneid] += 1
                fulloutdict[geneid][-2] = 1
            found2.append(geneid)
    
    for gene in genelist:
        if gene not in found2:
            fulloutdict[gene].append(0)
    print('Checking of alleles finished')
        
    print('genes with tissue annotation:', len(found))
    print('genes with allele causing phenotype:',len(found2))    

    return fulloutdict, tallyany, tallyphenotype

#Output files are made, listing the score associated with each gene based on the totality of available data
#If multiple terms have been provided, a gene will score 1 if associated with/causative of ANY OF THEM
def make_outputs(outfolder, outfoldertally, num_annotations, genelist, fulloutdict, tallyany, tallyphenotype):
    with open(f'{outfolder}/FullAnalysis.csv', 'w+') as fullanalysisfile:
        OUTfull = csv.writer(fullanalysisfile)

        OUTfull.writerow(fulloutdict['Header'])

        for gene in genelist:
            OUTfull.writerow(fulloutdict[gene])
    
    with open(f'{outfoldertally}/FlyBase_tally.csv', 'w+') as tallyoutfile:
        OUTtally = csv.writer(tallyoutfile)

        header = ['Gene ID', 'FlyBase - Score for ANY association with defined tissues', 'Flybase - Score for KNOWN PHENOTYPE in defined tissues']
        OUTtally.writerow(header)

        for gene in genelist:
            outrow = [gene, (tallyany[gene]/num_annotations), (tallyphenotype[gene]/num_annotations)]
            OUTtally.writerow(outrow)

def check_the_literature(input_fbids, target_annotations, outputfb, outputtally, alleles_file):

    annolist = target_annotations
    writable_annotations = {}
    
    #Controlled vocabulary terms are made machine readable
    for base_anno in target_annotations:
        i_anno = '' + base_anno

        if '/' in i_anno:
            i_anno = i_anno.replace('/', '_')
        if ' ' in i_anno:
            i_anno = i_anno.replace(' ', '_')

        writable_annotations[base_anno] = i_anno

    print(writable_annotations)
    
    #The genelist is populated
    genelist = populate_lists(input_fbids)

    num_annotations = len(annolist)

    tallyany = {}
    tallyphenotype = {}
    fulloutdict = {}

    fulloutdict['Header'] = ['Gene ID', ]
    
    #Output lists/dicts are initialised
    for gene in genelist:
        tallyany[gene] = 0
        tallyphenotype[gene] = 0
        fulloutdict[gene] = [gene,]
   
    #For each CV term, FlyBase is queried for associated genes
    #Gene tallies are then updated based on this data.
    for annotation in annolist:
        writable_anno = writable_annotations[annotation]
        SeleniumDriver = Web_Interface(annotation, outputfb, writable_anno)
        error_cases = SeleniumDriver.find_annotated_genes()

        fulloutdict, tallyany, tallyphenotype = check_lists(annotation, outputfb, fulloutdict, tallyany, tallyphenotype, genelist, alleles_file, error_cases, writable_anno)
        print(f'{annotation} finished')

    make_outputs(outputfb, outputtally, num_annotations, genelist, fulloutdict, tallyany, tallyphenotype)
