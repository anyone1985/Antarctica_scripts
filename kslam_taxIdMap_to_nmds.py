#! python3
# kslam_taxIdMap_to_nmds.py
# This script reads in the output(s) of kslam_xml_taxIdMap.py and formats a file that can be
# read in to VEGAN in R. Provided arguments are the file name(s), the output NMDS-ready tab-delimited file name,
# and the taxonomic level you want to summarise counts under.

import re, os, argparse

# Define regular expression for data extraction
taxline_regex = re.compile(r'\d+_+(\w+)$')
taxonomy = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

##### USER INPUT SECTION

usage = """%(prog)s reads in text file(s) resulting from the 
kslam_xml_taxIdMap.py script to produce a table that is amenable to NMDS.
This script may be provided with either the individual file names in the cmd line, or a text file
containing the file names on a new line per entry. Other necessary inputs are the location of a taxonomy
file (possibly produced by MetaPalette's generate_taxonomy_taxid.py script), the output file name,
and the taxonomy level to summarise abundances at (e.g., genus, class, etc.)
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("pos_input", nargs = "*", help="list k-SLAM tax id mapped file(s) in cmd line")
p.add_argument("-i", "--i", "-input", "--input", dest="txt_input",
                  help="optional ability to provide a text file which lists the k-SLAM tax id mapped file(s) instead of providing it in cmd line")
p.add_argument("-t", "--t", "-taxonomy", "--taxonomy", dest="taxonomy",
                  help="taxonomy file (possibly produced by MetaPalette's generate_taxonomy_taxid.py script)")
p.add_argument("-o", "--o", "-output", "--output", dest="output",
                  help="output file name")
p.add_argument("-l", "--l", "-level", "--level", dest="level", choices = taxonomy,
                  help="taxonomic level to summarise abundances under")
args = p.parse_args()

# Obtain data from arguments
outName = args.output
taxonomy_name = args.taxonomy
taxonomy_level = args.level

# Get file name(s) depending on whether they were provided in the cmd line or in a text file
try:
        if args.pos_input == [] and args.txt_input != None:
                txt_name = args.txt_input
                txt_contents = open(txt_name, 'r')
                kslam_files = []
                for line in txt_contents:
                        if line == '' or line == '\r\n' or line == '\n': continue
                        kslam_files.append(line.rstrip('\r\n'))
        elif args.pos_input != []:
                kslam_files = args.pos_input
        else:
                print('Something went wrong with finding the tax id mapped k-SLAM file(s).')
                print('It seems like you didn\'t provide any arguments that pointed towards the file(s).')
                print('Make sure to either positionally supply the individual file names in the cmd line, or provide a text file with these names using the -i tag.')
except:
        print('Something went wrong with finding the tax id mapped k-SLAM files.')
        if args.pos_input == []:
                print('It is likely that I couldn\'t find the text file you pointed at using the -i tag. Make sure the file name is correct, or provide the full directory to this file.')
        else:
                print('I am not sure what happened here. Double-check your arguments to make sure everything is correct.')
print(kslam_files)

# Find the index of the specified tax level
tax_index = taxonomy.index(taxonomy_level.lower())

# Load the taxonomy file and make a dictionary
taxonomy_file = open(taxonomy_name, 'r')
taxDict = {}
for line in taxonomy_file:
        if line.startswith('1_root\t'): continue
        line_split = line.rstrip('\n').split('\t')
        tax_format = line_split[2].replace('norank__1_root|', '')
        # Put entry in dictionary
        taxDict[line_split[1]] = tax_format

#### CORE LOOP
savedDicts = []
for filename in kslam_files:
        # Load the kslam file
        kslamFile = open(filename, 'r').read().split('\n')

        # Replace taxonomy ID with full NCBI tax format [and calculate normalisation number to make each abundance number add up to 100(%), and create a blank dictionary for later count adding]
        updated_idmap = []
        normalisation = 0
        countsDict = {}
        for line in kslamFile:
                if line.startswith('#') or line == '': continue
                line_split = line.split('\t')
                taxID = line_split[1]
                updated_line = [line_split[2]] + [taxDict[taxID]]
                updated_idmap.append(updated_line)
                # Normalisation
                normalisation += float(line_split[2])
                # Dictionary
                dicTax = taxDict[taxID]
                levels = list(reversed(dicTax.split('|')))
                # Get specified taxonomy level if present in this line
                if '|' + taxonomy_level + '__' in dicTax:
                        for level in levels:
                                if level.startswith(taxonomy_level + '_'):
                                        name = taxline_regex.findall(level)[0]
                                        countsDict.setdefault(name, 0)

                # Get the next lowest taxonomy level otherwise
                else:
                        breakAll = 0
                        # Look down for the next lowest level [only down to species level]
                        lowerRanks = taxonomy[tax_index+1:]
                        for rank in lowerRanks:
                                if '|' + rank + '__' in dicTax:
                                        for level in levels:
                                                if level.startswith(rank + '_'):
                                                        temp = level.split('_')
                                                        name = '_'.join(temp[3:])
                                                        countsDict.setdefault(name, 0)
                                                        breakAll = 1
                                                        break
                                if breakAll == 1:
                                        #breakAll = 0
                                        break
                        if breakAll == 0:
                                # Look up for the next highest level instead
                                higherRanks = list(reversed(taxonomy[:tax_index]))
                                for rank in higherRanks:
                                        if '|' + rank + '__' in dicTax:
                                                for level in levels:
                                                        if level.startswith(rank + '_'):
                                                                temp = level.split('_')
                                                                name = '_'.join(temp[3:])
                                                                countsDict.setdefault(name, 0)
                                                                breakAll = 1
                                                                break
                                        if breakAll == 1:
                                                break
                        if breakAll == 0:
                                # Just get the deepest rank if we get to this point [typically means the only rank is a norank]
                                for level in levels:
                                        temp = level.split('_')
                                        name = '_'.join(temp[3:])
                                        countsDict.setdefault(name, 0)
                                        break

        # Tally counts at specified level
        for line in updated_idmap:
                levels = list(reversed(line[1].split('|')))
                if '|' + taxonomy_level + '__' in line[1]:
                        for level in levels:
                                if level.startswith(taxonomy_level + '_'):
                                        name = '_'.join(level.split('_')[3:])
                                        countsDict[name] = countsDict[name] + float(float(line[0])/normalisation)
                                        break
                else:
                        breakAll = 0
                        # Look down for the next lowest level [only down to species level]
                        lowerRanks = taxonomy[tax_index+1:]
                        for rank in lowerRanks:
                                if '|' + rank + '__' in line[1]:
                                        for level in levels:
                                                if level.startswith(rank + '_'):
                                                        temp = level.split('_')
                                                        name = '_'.join(temp[3:])
                                                        countsDict[name] = countsDict[name] + float(float(line[0])/normalisation)
                                                        breakAll = 1
                                                        break
                                if breakAll == 1:
                                        break
                        if breakAll == 0:
                                # Look up for the next highest level instead
                                higherRanks = list(reversed(taxonomy[:tax_index]))
                                for rank in higherRanks:
                                        if '|' + rank + '__' in line[1]:
                                                for level in levels:
                                                        if level.startswith(rank + '_'):
                                                                temp = level.split('_')
                                                                name = '_'.join(temp[3:])
                                                                countsDict[name] = countsDict[name] + float(float(line[0])/normalisation)
                                                                breakAll = 1
                                                                break
                                        if breakAll == 1:
                                                break
                        if breakAll == 0:
                                # Just get the deepest rank if we get to this point [typically means the only rank is a norank]
                                for level in levels:
                                        temp = level.split('_')
                                        name = '_'.join(temp[3:])
                                        countsDict[name] = countsDict[name] + float(float(line[0])/normalisation)
                                        break

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
for i in range(len(kslam_files)):
    sampleName = kslam_files[i].rstrip('.kslam_taxidMapped')
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
