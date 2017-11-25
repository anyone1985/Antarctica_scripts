#! python3
# top_taxa_kslam_multicol.py
# This script was made quickly just to provide assistance with parsing
# columns in Microsoft Excel (or equivalent) produced by the nmds_convenience.py script.
# To use this you need to edit the underlying script to reflect the name of your metapalette-formatted
# taxonomy file, as well a the taxa of interest you are searching for.
# Use it by copying the column(s) produced by aforementioned script then press ENTER.
# It will count the proportion that this taxa contributes to the overall microbial abundance.

import os, pyperclip, re

taxonomy_name = 'mpal_tax_format.txt' ### EDIT

taxonomy_file = open(taxonomy_name, 'r')
taxDict = {}
for line in taxonomy_file:
        if line == '' or line == 'norank__1_root\n': continue
        line = line.rstrip('\n').replace('norank__1_root|', '').split('\t')[2]
        lastLevel = line.split('|')[-1]
        lastName = '_'.join(lastLevel.split('_')[3:])
        taxDict[lastName] = line

### HELPER SCRIPT
abund_regex = re.compile(r'([A-z0-9_;-]+?)\t(0\.\d{1,12})')
taxOfInterest = ['Bacteria']  ### EDIT
while True:
        try:
                print('___________________________________________')
                print('count!')
                but=input()
                table=pyperclip.paste()
                hits = abund_regex.findall(table)
                # Pull out taxa names
                summed = []
                overallAbundance = 0
                for hit in hits:
                        summed.extend([hit[0]])
                        overallAbundance += float(hit[1])
                summed = list(set(summed))
                # Sum taxa abundances
                sumDict = {}
                for i in range(len(summed)):
                        name = summed[i]
                        val = 0
                        for hit in hits:
                                if name in hit:
                                        val += float(hit[1])
                        sumDict[name] = val
                # Make list of the dictionary
                dicList = list(reversed(sorted(sumDict.items(), key = lambda x: x[1])))
                # Pull out relevant summed hits
                for tax in taxOfInterest:
                        summedCount = 0
                        for entry in dicList:
                                fullTax = taxDict[entry[0]]
                                if '_' + tax in fullTax:
                                        print(entry[0] + '\t' + str(entry[1]) + ' is ' + tax)
                                        #print(entry)
                                        summedCount += entry[1]
                        print(tax + ' make up ' + str(summedCount) + ' of total abundance (' + str(overallAbundance) +')')
                        print('Proportion = ' + str(summedCount/overallAbundance))
                        print('')
        except KeyboardInterrupt:
                quit()
        except:
                print('Uh oh. Something went wrong. Try again. If it keeps dying, figure out what\'s going on.')
