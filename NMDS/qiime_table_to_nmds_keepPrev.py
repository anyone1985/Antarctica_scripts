#! python3
# qime_table_to_nmds_keepPrev.py
# Simple script that will read in a QIIME table file (originating from the summarize_taxa_through_plots.py script or equivalent) and
# create a table that is amenable to NMDS. The arguments needed for this script the input and output file names.
# If you are choosing to look at genus level,  you should use the L6 file. Other levels should similarly use the appropriate L# file.
# Note that this script will find the lowest level possible in whatever file you feed it. That means a genus level NMDS may contain
# entries from higher levels if QIIME did not assign reads to genus level.
#
# This script is not meant to replace qime_table_to_nmds_keepPrev.py, it is just useful
# for use with the nmds_convenience.py script for making a sorted table

import re, os, argparse

regex = re.compile(r'(\w__[A-Za-z0-9-_\[\]\s]+?)(;|$)')

##### USER INPUT SECTION

usage = """%(prog)s reads in user specified QIIME taxa summary files
(typically obtained through output of summarize_taxa_through_plots.py) and
creates an output file that can be read by VEGAN in R.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "--i", "-input", "--input", dest="input",
                  help="QIIME table (L2-L6)")
p.add_argument("-o", "--o", "-output", "--output", dest="output",
                  help="output file name")
args = p.parse_args()

# Obtain data from arguments
fileName = args.input
outName = args.output

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

# Obtain counts from file
savedDicts = []
for i in range(1, len(qiimeFile[0].split('\t'))):
        sampleDict = {}
        for line in qiimeFile:
                # Skip header
                if line.startswith('#'):
                        continue
                # Handle unassigned
                if line.startswith('Unassigned'):
                        name = 'Unassigned'
                        line = line.split('\t')
                        sampleDict[name] = float(line[i])
                        continue
                # Remove 'other' labels
                #line = line.replace(';Other', '')
                sl = line.split('\t')
                # Get deepest node
                nodes = regex.findall(sl[0])
                name = nodes[-1][0]
                # Added: get other nodes up to class
                extra = ''
                for node in reversed(nodes):
                        if node[0] == name:
                                continue
                        if not node[0].startswith('p_'):
                                extra = node[0].replace('__', '_') + ';' + extra
                        else:
                                break
                name = extra + name
                if name not in sampleDict:
                        sampleDict[name] = float(sl[i])
                else:
                        sampleDict[name] += float(sl[i])
        # Add into the overall dictionary list
        savedDicts.append(sampleDict)

#### Format an output table from the saved dictionaries

# Pull out all entries found in the data
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
