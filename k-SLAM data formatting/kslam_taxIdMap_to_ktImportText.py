#! python3
# kslam_taxIdMap_to_ktImportText.py
# This script reads in the output(s) of kslam_xml_taxIdMap.py and formats file(s) that can be
# called by ktImportText to produce Krona plot(s). Provided arguments are the file name(s), the output NMDS-ready tab-delimited file,
# and the taxonomic level you want to summarise counts under.

import re, os, argparse

# Define regular expression for data extraction
taxline_regex = re.compile(r'\d+_+(\w+)$')

##### USER INPUT SECTION

usage = """%(prog)s reads in text file(s) resulting from the 
kslam_xml_taxIdMap.py script to produce file(s) that can be called by ktImportText.
This script may be provided with either individual file name(s) in the cmd line, or a text file
containing the file name(s) on a new line per entry. The other input requires is the location of a full taxonomy
file (possibly produced by MetaPalette's generate_taxonomy_taxid.py script). Output file names will automatically
be made by adding '_kronaReady' before the file's original extension.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("pos_input", nargs = "*", help="list k-SLAM tax id mapped file(s) in cmd line")
p.add_argument("-i", "--i", "-input", "--input", dest="txt_input",
                  help="optional ability to provide a text file which lists the k-SLAM tax id mapped file(s) instead of providing it in cmd line")
p.add_argument("-t", "--t", "-taxonomy", "--taxonomy", dest="taxonomy",
                  help="full taxonomy file (possibly produced by MetaPalette's generate_taxonomy_taxid.py script)")
args = p.parse_args()

# Obtain data from arguments
taxonomy_name = args.taxonomy

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

# Load the taxonomy file and make a dictionary
taxonomy_file = open(taxonomy_name, 'r')
taxDict = {}
for line in taxonomy_file:
        line_split = line.rstrip('\n').split('\t')
        tax_format = line_split[2].lstrip('norank__1_root|')
        # Put entry in dictionary
        if line_split[1] not in taxDict:
                taxDict[line_split[1]] = tax_format
        else:
                if len(tax_format) > len(taxDict[line_split[1]]):
                        taxDict[line_split[1]] = tax_format

#### CORE LOOP
end_print = []
for filename in kslam_files:
        kslamFile = open(filename, 'r').read().split('\n')
        # Replace taxonomy ID with full NCBI tax format
        updated_idmap = []
        for line in kslamFile:
                if line.startswith('#') or line == '': continue
                line_split = line.split('\t')
                taxID = line_split[1]
                updated_line = [line_split[2]] + [taxDict[taxID]]
                updated_idmap.append(updated_line)

        # Break taxonomy lines into krona-ready format
        outputList = []
        for line in updated_idmap:
                taxLevels = line[1].split('|')
                updatedTaxLevels = []
                for level in taxLevels:
                        name = taxline_regex.findall(level)[0]
                        updatedTaxLevels.append(name)
                updatedTaxLevels = '\t'.join(updatedTaxLevels)
                outputList.append(str(float(line[0])*100) + '\t' + updatedTaxLevels)
                
        # Output name
        suffix = filename.split('.')[-1]
        prefix = '.'.join(filename.split('.')[0:-1])
        outName = prefix + '_kronaReady.' + suffix
        end_print.append(outName)
        
        # Create output file
        output = '\n'.join(outputList)
        outFile = open(outName, 'w')
        outFile.write(output)
        outFile.close()


#########
print('Copy paste the below text if you want to use this for the ktImportText program call:')
print(' '.join(end_print))
