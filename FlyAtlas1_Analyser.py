import csv
import os

#Creates a dictionary associating each sequence in the microarray with specific gene(s)
def populate_annotations(annotationfile):

    with open(annotationfile) as annofile:
        INanno = csv.reader(annofile)
        next(INanno)

        annotation_dict = {}

        for line in INanno:
            if line[0][0] != '#' and line[0] != 'Probe Set ID':
                annotation = line[0] 
                flybaseid = line[24]

                annotation_dict[annotation] = flybaseid
    
    return(annotation_dict)

#Gathers tissues of interest from designated text files, and associated them with the relevant tissue identifies from FlyAtlas1
def build_tiss_lists(unconvertedlist):
    fulllistadult = []
    fulllistlarval = []
    tcdictadult = {}
    tcdictlarval = {}

    #No-one's favorite, but a big list containing all tissues and their exact designation in the FlyAtlas1 dataset is necessary due to inconsistent naming conventions for different tissues
    tissueconversions = [('Brain/CNS', 'Brain', 'Lcns'),('Head', 'head', 'n/a'), ('Crop', 'crop', 'n/a'), ('Midgut', 'Midgut', 'l_mid'), ('Hindgut', 'Hindgut', 'l_hind'), ('Malpighian tubule','tubule', 'l_tub'), ('Ovary', 'Ovary', 'n/a'),
    ('Testis', 'Testis', 'n/a'), ('Accessory gland', 'Acc', 'n/a'), ('Fat Body', 'AFB', 'l_fat'), ('Thoracicoabdominal ganglion', 'tag', 'n/a'), ('Carcass', 'car', 'Lcar'), ('Salivary gland', 'Sg', 'l_sg'), 
    ('Virgin spermatheca', 'SptV', 'n/a'), ('Mated spermatheca', 'SptM', 'n/a'), ('Eye', 'eye', 'n/a'), ('Heart', 'heart', 'n/a'), ('Trachea', 'n/a', 'trachea'), ('S2', 'n/a', 'n/a')]

    for pair in tissueconversions:
        tcdictadult[pair[0]] = pair[1]
        tcdictlarval[pair[0]] = pair[2]
    
    for item in unconvertedlist:
        if tcdictadult[item] != 'n/a':
            fulllistadult.append(tcdictadult[item])
        if tcdictlarval[item] != 'n/a':
            fulllistlarval.append(tcdictlarval[item])

    return fulllistadult, fulllistlarval

#The centrepiece function, which takes raw abundances from the FlyAtlas1, uses them to generate enrichment values based on <AVERAGES>, see below, and then reports genes expressed in tissues of
#in adult flies, larval flies, and both.
#This can run in three modes - Default mode, in which enrichments are simply extracted from the FlyAtlas1 dataset (tissue specific/whole fly)
#                            - Relative mode in which enrichments are based on adult abundances/average of adult tissues and larval abundances/average of larval tissues
#                            - All mode, in which both adult and larval tissue abundances are divided by the average reads across all tissues

def get_expression(annotation_dict, enr_thresh, abn_thresh, mode, adult_toi, larval_toi, presenttolerance, output, fa1data, permissive_adult, permissive_larval):

    with open(fa1data) as atlas1file:
        atlas1IN = csv.reader(atlas1file, delimiter = '\t')

        linecount = 0
        tissues = {}
        readdict = {}

        enrich_adult = []
        enrich_larval = []
        enrich_all = []

        abun_adult = []
        abun_larval = []
        abun_alsep = []

        for line in atlas1IN:
            adult_enrich = 0
            larval_enrich = 0
            all_enrich = 0

            adultreads = []
            larvalreads = []
            allreads = []

            toi_abundances_adult = []
            other_abundances_adult = []

            toi_abundances_larval = []
            other_abundances_larval = []

            probeid = line[0]

            adult_enrichments = []
            larval_enrichments = []

            updatedavgs_adult = []
            updatedavgs_larval = []
            updatedavgs_all = []

            #Only grabs the "Mean" column from each tissue in the .tsv
            for index in range(2, 130, 5):
                
                #builds a dictionary defining where each tissue is within the FlyAtlas1 dataset
                if linecount == 0:

                    tissuename = line[index]
                    tissuenamefix = tissuename[:-4]
                    tissuenamefix = tissuenamefix.strip(' ')
                    tissues[index] = tissuenamefix
          

                else:
                    #Skips all rows without a valid probe ID
                    if probeid != '':
                        #Gives each probe the appropriate FbID annotation and grabs useful digits from the row
                        annotated = annotation_dict[probeid]
                        presence = index + 2
                        ratio = index + 3
                        tissuetype = tissues[index]

                        #Genes are only considered for further analysis if they meet or exceed the user-defined present call tolerance in the tissues of interest
                        if int(line[presence]) >= presenttolerance:
                            
                            #Reads from tissues on the Tissue of Interest lists are held in readdict, for ease of manipulation later. 
                            if tissuetype in adult_toi:
                                readdict[tissuetype] = float(line[index])
                                adultreads.append(float(line[index]))
                                allreads.append(float(line[index]))
                                
                                toi_abundances_adult.append(float(line[index]))

                                adult_enrichments.append(float(line[ratio]))

                            elif tissuetype in larval_toi:
                                readdict[tissuetype] = float(line[index])
                                larvalreads.append(float(line[index]))
                                allreads.append(float(line[index]))
                                
                                toi_abundances_larval.append(float(line[index]))

                                larval_enrichments.append(float(line[ratio]))

                            #Background abundances (ie abundances from non-target tissues) are caught to check for increased abundance in tissues of interest
                            else:
                                #We need to catch the tissues that are larval to split these off from adult tissues 
                                if tissuetype[0] == 'l' or tissuetype[0] == 'L' or tissuetype == 'trachea ':
                                    larvalreads.append(float(line[index]))
                                    allreads.append(float(line[index]))
                                    
                                    if tissuetype not in permissive_larval:
                                        other_abundances_larval.append(float(line[index]))
                                
                                else:
                                    adultreads.append(float(line[index]))
                                    allreads.append(float(line[index]))
                                
                                    if tissuetype not in permissive_adult:
                                        other_abundances_adult.append(float(line[index]))
                        
                        #If the present call falls below the defined threshold, the gene is assumed to not be expressed in a tissue.
                        #This isn't strictly accurate but fulfils the necessary purpose 
                        # -> genes which aren't consistently expressed in target tissues are not pulled out as genes of interest
                        else:
                            readdict[tissuetype] = 0

            if probeid != '' and linecount != 0:
                #Finds the average expression of each gene in all adult tissues...
                try:
                    adultavg = (sum(adultreads)/len(adultreads))
                except:
                    adultavg = 'N/A'
                
                #All larval tissues...
                try:
                    larvalavg = (sum(larvalreads)/len(larvalreads))
                except:
                    larvalavg = 'N/A'
                
                #And ALL tissues.
                try:
                    allavg = (sum(allreads)/len(allreads))
                except:
                    allavg = 'N/A'
                
                #Runs default mode - enrichment ratios are simply taken from the FlyAtlas1 .tsv and compared against the defined threshold
                if mode == 'D' or mode == 'd':
                    #This corrects for cases where a failure to meet the presentcount threshold results in background lists which are completely empty
                    while len(adult_enrichments) < len(adult_toi):
                        adult_enrichments.append(0)
                    while len(larval_enrichments) < len(larval_toi):
                        larval_enrichments.append(0)

                    try:
                        if min(adult_enrichments) >= enr_thresh:
                            enrich_adult.append(annotated)
                    except:
                        pass
                    
                    try:
                        if min(larval_enrichments) >= enr_thresh:
                            enrich_larval.append(annotated)
                    except:
                        pass

                    if annotated in enrich_adult and annotated in enrich_larval:
                        enrich_all.append(annotated)

                else:
                    #Runs ALL mode - the average gene abundance in ALL tissues is used to generate enrichment values to compare against threshold
                    if mode == 'A' or mode == 'a':
                        divisor_adult = allavg
                        divisor_larva = allavg
                    
                    #Runs RELATIVE mode - the relative average abundance in adult and larval tissue is used to generate enrichment values to compare against threshold
                    elif mode == 'R' or mode == 'r':
                        divisor_adult = adultavg
                        divisor_larva = larvalavg

                    for tissue in adult_toi:
                        try:
                            newratio = readdict[tissue]/divisor_adult
                            updatedavgs_adult.append(newratio)
                            updatedavgs_all.append(newratio)
                            if newratio >= enr_thresh:
                                adult_enrich += 1
                                all_enrich += 1
                        except:
                            updatedavgs_adult.append(0)
                            updatedavgs_all.append(0)

                    for tissue in larval_toi:
                        try:
                            newratio = readdict[tissue]/divisor_larva
                            updatedavgs_larval.append(newratio)
                            updatedavgs_all.append(newratio)
                            if newratio >= enr_thresh:
                                larval_enrich += 1
                                all_enrich += 1
                        except:
                            updatedavgs_larval.append(0)
                            updatedavgs_all.append(0)

                #This corrects for cases where a failure to meet the presentcount threshold results in background lists which are completely empty
                if other_abundances_adult == []:
                    other_abundances_adult.append(0)
                if other_abundances_larval == []:
                    other_abundances_larval.append(0)
                
                #Again correcting for cases where presentcount thresholds result in shortened lists, in this case for tissues of interest 
                while len(toi_abundances_adult) < len(adult_toi):
                    toi_abundances_adult.append(0)
                while len(toi_abundances_larval) < len(larval_toi):
                    toi_abundances_larval.append(0)

                #Here, we find our genes of interest - genes which are more abundant in tissue of interest lists vs background lists are gathered for both fly types individually and for genes at increased abundance in both
                #Similarly, we find genes which are expressed >= the enrichment threshold in all tissues of interest in each fly type individually and in both types.
                if linecount > 0 and probeid != '':
                    try:
                        if min(toi_abundances_adult) > abn_thresh * max(other_abundances_adult):
                            abun_adult.append(annotated)

                            try:
                                if min(toi_abundances_larval) > abn_thresh * max(other_abundances_larval):
                                    abun_alsep.append(annotated)
                            except:
                                pass
                    except:
                        pass

                    try:
                        if min(toi_abundances_larval) > abn_thresh * max(other_abundances_larval):
                            abun_larval.append(annotated)
                    except:
                        pass
                    
                    if mode != 'D' and mode != 'd':
                        if adult_enrich == len(adult_toi):
                            enrich_adult.append(annotated)
                        if larval_enrich == len(larval_toi):
                            enrich_larval.append(annotated)
                        if all_enrich == (len(adult_toi) + len(larval_toi)):
                            enrich_all.append(annotated)
                    
            linecount += 1

    #Here, we print seperate output ABUNDANCE and ENRICHMENT gene of interest files for each fly type, and a third set for "ALL stages"
    for stage in ['ADULTS', 'LARVAL', 'ALL']:

        printed_e = []
        printed_a = []

        if stage == 'ADULTS':
            enrichlist = enrich_adult
            abunlist = abun_adult

        elif stage == 'LARVAL':
            enrichlist = enrich_larval
            abunlist = abun_larval

        elif stage == 'ALL':
            enrichlist = enrich_all
            abunlist = abun_alsep

        with open(f'{output}/{stage}_ENRICHED.txt', 'w+') as enrichfile, open(f'{output}/{stage}_ABUNDANT.txt', 'w+') as abunfile:
            
            #the If statements below just account for cases where probes could not be adequately annotated earlier with an FbID.
            for item in enrichlist:
                if item != '---' and item not in printed_e:
                    enrichfile.write(f'{item}\n')
                    printed_e.append(item)

            for item in abunlist:
                if item != '---' and item not in printed_a:
                    abunfile.write(f'{item}\n')
                    printed_a.append(item)

#Converts the float values for thresholds to strings, replacing "."s with "-"s for filenames
def make_string(floating_point, identifier):
    string_vers = str(floating_point)
    string_vers = string_vers.replace('.', '-')

    final_string = identifier + string_vers

    return final_string

#Updates the output folder to reflect enrichment thresholds, then creates the folder
#def update_out_name(folder_name, enr_thresh, abn_thresh):
    #enrthresh_string = make_string(enr_thresh,'enr')
    #abnthresh_string = make_string(abn_thresh, 'abn')

    #new_out_name = folder_name + '_' + enrthresh_string + '_' + abnthresh_string

    #os.mkdir(f'{new_out_name}')

    #return new_out_name

def analyse_FA1(annotationfile, targettissues, permissive, outfolder, enrstringency, abnstringency, mode, presentthresh, fa1file):
    
    annodict = populate_annotations(annotationfile)

    toiadult, toilarval = build_tiss_lists(targettissues)

    permissive_adult = []
    permissive_larval = []

    if permissive != [] and permissive != ['']:
        permissive_adult, permissive_larval = build_tiss_lists(permissive)

    #new_output_folder = update_out_name(outfolder, enrstringency, abnstringency)

    get_expression(annodict, enrstringency, abnstringency, mode, toiadult, toilarval, presentthresh, outfolder, fa1file, permissive_adult, permissive_larval)
