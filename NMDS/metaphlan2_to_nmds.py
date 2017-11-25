#! python3
# metaphlan2_to_nmds.py
# This script reads in the output(s) of metaphlan2 and formats a file that can be
# read in to VEGAN in R. Provided arguments are the file name(s), the output NMDS-ready tab-delimited file name,
# and the taxonomic level you want to summarise counts under.

import os, argparse

# Define regular expression for data extraction
taxonomy = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus']
taxAbbrev = {'superkingdom': 'k', 'phylum': 'p', 'class': 'c', 'order': 'o', 'family': 'f', 'genus': 'g'}

##### USER INPUT SECTION

usage = """%(prog)s reads in the text file(s) resulting from MetaPhlAn2
and returns taxa counts at the specified taxonomic level. Taxonomic levels
are restricted since MetaPhlAn2 is not guaranteed to assign below genus level.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("pos_input", nargs = "*", help="list MetaPhlAn2 results file(s) in cmd line")
p.add_argument("-i", "-input", dest="txt_input",
                  help="optional ability to provide a text file which lists the MetaPhlAn2 results file(s) instead of providing it in cmd line")
p.add_argument("-o", "-output", dest="output",
                  help="output file name")
p.add_argument("-l", "-level", dest="level", choices = taxonomy,
                  help="taxonomic level to summarise abundances under")

args = p.parse_args()

# Obtain data from arguments
outName = args.output
taxonomy_level = args.level

# Get file name(s) depending on whether they were provided in the cmd line or in a text file
if args.txt_input != None:
        txt_name = args.txt_input
        txt_contents = open(txt_name, 'r')
        metaphlan2_files = []
        for line in txt_contents:
                if line == '' or line == '\r\n' or line == '\n': continue
                metaphlan2_files.append(line.rstrip('\r\n'))
        # Check if files exist
        for file in metaphlan2_files:
                if not os.path.isfile(file):
                        print('Can\'t find ' + file + '. Make sure the name is correct, or you specify the path to this file if it is not in the current directory')
                        quit()
elif args.pos_input != []:
        metaphlan2_files = args.pos_input
        # Check if files exist
        for file in metaphlan2_files:
                if not os.path.isfile(file):
                        print('Can\'t find ' + file + '. Make sure the name is correct, or you specify the path to this file if it is not in the current directory')
                        quit()
elif args.pos_input == [] and args.txt_input == None:
        print('Something went wrong with finding the MetaPhlAn2 file(s).')
        print('It seems like you didn\'t provide any arguments that pointed towards the file(s).')
        print('Make sure to either positionally supply the individual file names in the cmd line, OR provide a text file with these names using the -i tag.')

print(metaphlan2_files)

taxonomy_level = taxAbbrev[taxonomy_level]

#### CORE LOOP
savedDicts = []
for filename in metaphlan2_files:
        # Load the metaphlan2 file
        metaphlan2File = open(filename, 'r')
        # Tally counts at specified level
        countsDict = {}
        for line in metaphlan2File:
                if '|' + taxonomy_level + '__' not in line:
                        continue
                line=line.rstrip('\n').split('\t')
                taxLevels = line[0].split('|')
                if not taxLevels[-1].startswith(taxonomy_level + '__'):
                        continue
                # Add taxa to the dictionary
                name = taxLevels[-1].split(taxonomy_level + '__')[1]
                countsDict.setdefault(name, 0)
                countsDict[name] = countsDict[name] + float(float(line[1])/100) # Make it so it sums to 1.0, not 100.0
        # Save dictionary as a list for later sorting
        savedDicts.append(countsDict)

#### Format an output table from the saved dictionaries

# Pull out all entries found in the data at the specified taxonomic level
taxonomyList = []
for dictionary in savedDicts:
    for entry in dictionary.items():
        taxonomyList.append(entry[0])
taxonomyList = sorted(list(set(taxonomyList)))

# Tabulate results
header = taxonomyList[:]
header.insert(0, '')
table = []
table.append(header)
for i in range(len(metaphlan2_files)):
    sampleName = metaphlan2_files[i].split('_')[1].upper()
    row = [sampleName]
    for entry in taxonomyList:
        if entry in savedDicts[i]:
            row.append(str(savedDicts[i][entry]))
        else:
            row.append('0')
    table.append(row)

# Output file
for x in range(len(table)):
    table[x] = '\t'.join(table[x])
output = '\n'.join(table)
outFile = open(outName, 'w')
outFile.write(output)
outFile.close()

print('Finished')
