#! python3
# kslam_xml_to_alphaCounts
# Script that will read in a kslam output .xml file and extract the actual mapped reads count
# and create a file that is amenable to further scripts for reformatting (intended for kslam_alphaCounts_to_simpson)

import re, os, argparse

# Define regular expression for data extraction
kslam_regex = re.compile(r'<taxon>\n.*?<abundance numReads="(\d{1,10})">(\d\.\d{1,10})<\/abundance>\n.*?<taxonomyID>(\d{1,10})<\/taxonomyID>\n.*?\n.*?<name>(.*?)<\/name>', re.DOTALL)

##### USER INPUT SECTION

usage = """%(prog)s reads in a k-SLAM xml output file and extracts the
actual mapped reads count, creating a file amenable to further scripts for reformatting
(intended for kslam_alphaCounts_to_simpson.py to follow this for alpha diversity calculation)

(Note that the output file will be called the input name + '_alphaCounts' and will be generated in the same directory as the input file)
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("xml", help="k-SLAM xml file")

args = p.parse_args()

# Obtain data from arguments
kslam_xml = args.xml
outName = kslam_xml + '_alphaCounts'

# Load the xml file and parse its contents
kslamFile = open(kslam_xml, 'r').read()

# Pull results with regex
regexHits = kslam_regex.findall(kslamFile)

# Parse results
outputList = ['#speciesName\t#taxID\t#read_counts']
for hit in regexHits:
        if hit[3] == '' and hit[2] == '31508':                  # Need to add this here since k-SLAM did not add a name in the xml where there was a hit with tax ID of 31508, and the taxonomy file produced by the MetaPalette script lacked this entry (has been collapsed into tax ID 1962501, likely reason for this malfunction)
                continue
        else:
                line = hit[3] + '\t' + hit[2] + '\t' + hit[0]
        outputList.append(line)

# Output file
output = '\n'.join(outputList)
outFile = open(outName, 'w')
outFile.write(output)
outFile.close()
