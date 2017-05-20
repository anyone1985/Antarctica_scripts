#! python3
# qiime_to_ktImportText.py
# Simple script that will read in a QIIME table file (originating from the summarize_taxa_through_plots.py script or equivalent) and
# create a text file that is amenable to ktImportText. The only argument needed for this script is the input file name.
# If you are choosing to look at genus level,  you should use the L6 file. Other levels should similarly use the appropriate L# file.

import re, os, argparse

# Define regular expression for data extraction
taxline_regex = re.compile(r'__(\w+)$')

##### USER INPUT SECTION

usage = """%(prog)s reads in user specified QIIME taxa summary files
(typically obtained through output of summarize_taxa_through_plots.py) and
creates an output text file that can be fed into ktImportText to produce Krona plots.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "--i", "-input", "--input", dest="input",
                  help="QIIME table sorted L6 file")
args = p.parse_args()

# Obtain data from arguments
qiimeName = args.input
dirName = os.getcwd() + '\\script_out\\'

# Create output dir if it does not already exist
if not os.path.isdir(dirName):
        os.mkdir(dirName)

# Load the qiime file
qiimeFile = open(qiimeName, 'r').read().split('\n')
if qiimeFile[0].startswith('# Constructed from biom file'):
        del qiimeFile[0]
if qiimeFile[-1] == '':
        del qiimeFile[-1]

# Break file into columns
for i in range(len(qiimeFile)):
        qiimeFile[i] = qiimeFile[i].split('\t')

# Break taxonomy lines into krona-ready format
outputList = []
end_print = ''
for i in range(1, len(qiimeFile[0])):
        for row in qiimeFile:
                if row[0] == '#OTU ID': continue        # Skip first line
                if float(row[i])*100 == 0: continue     # Skip null groups
                if row[0].startswith('Unassigned'):     # Skip name processing for unassigned group
                        outputList.append(str(float(row[i])*100) + '\tUnassigned')
                        continue
                taxLevels = row[0].split(';')
                updatedTaxLevels = []
                # Process the taxonomy level names
                for level in taxLevels:
                        if level == 'Other':
                                continue
                        name = taxline_regex.findall(level)
                        if len(name) == 0:
                                continue
                        updatedTaxLevels.append(name[0])
                updatedTaxLevels = '\t'.join(updatedTaxLevels)
                # Output
                outputList.append(str(float(row[i])*100) + '\t' + updatedTaxLevels)
                
        # Output
        output = '\n'.join(outputList)
        #outFile = open(outName + '_' + qiimeFile[0][i] + '.txt', 'w')
        outFile = open(dirName + qiimeFile[0][i] + '.txt', 'w')
        outFile.write(output)
        outFile.close()
        outputList = []
        end_print += qiimeFile[0][i] + '.txt  '

print('Copy paste the below text if you want to use this for the ktImportText program call:')
print(end_print)
