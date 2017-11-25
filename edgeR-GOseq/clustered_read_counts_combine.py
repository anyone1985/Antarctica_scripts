#! python3
# clustered_read_counts_combine.py
# This program reads in a group of clustered read counts files originating from the cluster_read_counts.py script and converts this into a table that is readable by edgeR
# As with the previous script, this could be memory intensive depending on file sizes

# Import necessary packages
import os, argparse

#### USER INPUT SECTION
usage = """%(prog)s reads in a group of clustered read counts files
 originating from the cluster_read_counts.py script and converts this into a table that is readable by edgeR.
 The header of this file will have columns labelled according to each filename up to the first '.'
 (e.g., 'SITE1.clstrd.counts' will have its column labelled 'SITE1').
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-d", "-directory", dest="directory",
                  help="Enter the directory where the clustered htseq read counts files can be located. Enter '.' to signify the current directory.")
p.add_argument("-e", "-extension", dest="extension",
                  help="Enter the extension of the clustered htseq read counts files (e.g., '.clstrd.counts') to distinguish which files this script should read. This should be unique to the files originating from cluster_read_counts.py.")
p.add_argument("-o", "--o", "--output", "-output", dest="output",
             help="output file name")
args = p.parse_args()

# Obtain data from arguments
directory = args.directory
extension = args.extension
outName = args.output

# Obtain data from htseq output files
fileNames = []
for file in os.listdir(directory):
        if file.endswith(extension):
                fileNames.append(os.path.join(directory, file))

#### CORE LOOP
ongoingCount = 0
for i in range(len(fileNames)):
        if i == 0:
                htseqFile = open(fileNames[i], 'r').read()
                combinedList = htseqFile.split('\n')
                while combinedList[-1] == '' or combinedList[-1].replace('\t', '') == '':               # Deal with blank lines if they exist. They shouldn't by default, but there's a few ways these could be introduced unintentionally
                        del combinedList[-1]
                for i in range(len(combinedList)):
                        combinedList[i] = combinedList[i].split('\t')
                print('Finished file (' + fileNames[ongoingCount] + ')')
                htseqFile = ''
        else:
                htseqFile = open(fileNames[i], 'r').read().split('\n')
                while htseqFile[-1] == '' or htseqFile[-1].replace('\t', '') == '':
                        del htseqFile[-1]
                htseqCounts = []
                for row in htseqFile:
                        htseqCounts.append(row.split('\t')[1])
                for i in range(len(htseqFile)):
                        combinedList[i].append(htseqCounts[i])
                print('Finished file (' + fileNames[ongoingCount] + ')')
                htseqFile = ''
        ongoingCount += 1

# Output file
sampleNames = ['gene_id']
for name in fileNames:
        sampleNames.append(name.split('.')[0])
header = '\t'.join(sampleNames) + '\n'
for i in range(len(combinedList)):
        combinedList[i] = '\t'.join(combinedList[i])
output = '\n'.join(combinedList)
output = header + output
outFile = open(outName, 'w')
outFile.write(output)
outFile.close()
print('All done')
