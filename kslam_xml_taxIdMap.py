#! python3
# kslam_xml_taxIdMap.py
# This script will read in the kslam xml output file and create a tab-delimited text
# file that can be used for further downstream scripts that can create a text file
# amenable to ktImportText Krona plot generation, or be used for NMDS plot generation using VEGAN.

import re, os, argparse

# Define regular expression for data extraction
kslam_regex = re.compile(r'<taxon>\n.*?<abundance numReads="(\d{1,10})">(\d\.\d{1,10})<\/abundance>\n.*?<taxonomyID>(\d{1,10})<\/taxonomyID>\n.*?\n.*?<name>(.*?)<\/name>', re.DOTALL)

##### USER INPUT SECTION

usage = """%(prog)s reads in a k-SLAM xml file and formats
something that can be called by ktImportTaxonomy to produce a Krona plot
or make an NMDS plot using VEGAN
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "--i", "-input", "--input", dest="xml",
                  help="k-SLAM xml file")
p.add_argument("-o", "--o", "-output", "--output", dest="output",
                  help="output file name")
args = p.parse_args()

# Obtain data from arguments
kslam_xml = args.xml
outName = args.output

# Load the xml file and parse its contents
kslamFile = open(kslam_xml, 'r').read()

# Pull results with regex
regexHits = kslam_regex.findall(kslamFile)

# Parse results
outputList = ['#speciesName\t#taxID\t#rel_abundance']
for hit in regexHits:
        if hit[3] == '' and hit[2] == '31508':          # Needed to add this here since k-SLAM did not add a name in the xml where there was a hit with tax ID of 31508, and the taxonomy file produced by the MetaPalette script lacked this entry (has been collapsed into tax ID 1962501, likely reason for this malfunction). This taxa was predicted at very low abundance in this study's results, so removing it was considered the best option.
                continue
        else:
                line = hit[3] + '\t' + hit[2] + '\t' + hit[1]
        outputList.append(line)

# Output file
output = '\n'.join(outputList)
outFile = open(outName, 'w')
outFile.write(output)
outFile.close()
