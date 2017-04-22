#! python3
# qime_table_to_nmds.py
# Simple script that will read in a QIIME table file (originating from the taxa_summary_plots script or equivalent) and
# create a table that is amenable to NMDS. The arguments needed for this script is the level that this script is looking at
# (e.g., 'c' means class, 'g' means genus) as well as the input and output file names. If you are choosing to look at genus or 'g' level
# you should use the L6 file. Other levels should similarly use the appropriate L# file.

import re, os, argparse

##### USER INPUT SECTION

usage = """%(prog)s reads in user specified k-SLAM files that have been
run through the qiime_xml_taxIdMap.py script to produce a table that
is amenable to NMDS. After the -i tag, list all the files with spaces in between their names
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "--i", "-input", "--input", dest="input",
                  help="QIIME table (L2-L6)")
p.add_argument("-o", "--o", "-output", "--output", dest="output",
                  help="output file name")
p.add_argument("-l", "--l", "-level", "--level", dest="level",
                  help="taxonomy level letter (e.g., enter 'c' without quotation marks for class, 'g' for genus)")
args = p.parse_args()

# Obtain data from arguments
fileName = args.input
outName = args.output
taxonomy_level = args.level

#### CORE LOOP
savedDicts = []

# Load the qiime file
qiimeFile = open(fileName, 'r').read().split('\n')
if qiimeFile[0].startswith('# Constructed from biom file'):     # Remove useless row from start of file
        del qiimeFile[0]
while qiimeFile[-1] == '':                                      # Remove blank rows from end of file
        del qiimeFile[-1]
# Pull out qiime header with sample IDs for later use
qiimeFileHeader = qiimeFile[0].split('\t')
if qiimeFileHeader[0] == '#OTU ID':
        del qiimeFileHeader[0]

# Tally counts at specified level
savedDicts = []
for i in range(1, len(qiimeFile[0].split('\t'))):
        sampleDict = {}
        for line in qiimeFile:
                # Handle unassigned
                if 'Unassigned' in line:
                        name = 'Unassigned'
                        line = line.split('\t')
                        sampleDict[name] = float(line[i])
                # Check if this row contains the taxonomy level
                if ';' + taxonomy_level + '__' not in line:
                        continue
                line = line.split('\t')
                taxLevels = line[0].split(';')
                # Add to sample dictionary
                for level in taxLevels:
                        if level.startswith(taxonomy_level + '_') and level != taxonomy_level + '__':
                                name = level[3:]
                                sampleDict[name] = float(line[i])
        # Add into the overall dictionary list
        savedDicts.append(sampleDict)

#### Format an output table from the saved dictionaries

# Pull out all entries found in the data at the specified taxonomic level
taxonomyList = []
for dictionary in savedDicts:
    for entry in dictionary.items():
        taxonomyList.append(entry[0])
taxonomyList = sorted(list(set(taxonomyList)))

# Tabulate results
header = taxonomyList[:]
header.insert(0, 'OTU ID')
table = []
table.append(header)
for i in range(0, len(qiimeFileHeader)):
    sampleName = qiimeFileHeader[i]
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
