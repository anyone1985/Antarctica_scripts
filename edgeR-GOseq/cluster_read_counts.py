#! python3
# cluster_read_counts.py
# This program reads in a cdhit .clstr file and htseq output file and clusters reads from redundant mappings into a single count corresponding to the representative sequence from each cluster
# This script is intended as a precursor to the clustered_read_counts_combine.py script which will take individually clustered files and create a table that edgeR can read.
# This script is likely to be a bit heavy on memory depending on file size. 

# Import necessary packages
import re, os, time, argparse

# Set up global expressions for later use
outputText = []
numReg = re.compile(r'\d{0,10}')
transIdReg = re.compile(r'>(\w{7,13})\.\.\.\s(\*)?')

#### USER INPUT SECTION
usage = """%(prog)s reads in cdhit .clstr file and htseq output file
and clusters reads from redundant mappings into a single count
corresponding to the representative sequence from each cluster
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-c", "--c", "-clstr", "--clstr", dest="clstr",
                  help=".clstr file")
p.add_argument("-ht", "--ht", "--htseq", "-htseq", dest="htseq",
                   help="htseq file")
p.add_argument("-o", "--o", "--output", "-output", dest="output",
             help="output file name")

args = p.parse_args()

# Obtain data from arguments
fileName = args.clstr
fileName2 = args.htseq
outName = args.output

# Load the clstr file and parse its contents
clstrFile = open(fileName, 'r').read()
chunkedClstr = clstrFile.split('\n>')
print('Finished chunking clustr file')

# Load the htseq output file and parse its contents
htseqFile = open(fileName2, 'r').read()
htseqList = htseqFile.split('\n')
htseqFile = ''  # Dump from memory
if htseqList[-1] == '':
        del htseqList[-1]
htseqDict = {}
for row in htseqList:
        temp = row.split('\t')
        htseqDict[temp[0]] = temp[1]
print('Made the dict')

#### CORE LOOP
start_time=time.time()
ongoingCount = 0
for cluster in chunkedClstr:
# Pull out sequence IDs from cluster and note the represenative sequence
        idsRegPull = transIdReg.findall(cluster)
        transIds = []
        for i in range(len(idsRegPull)):
                transIds.append(idsRegPull[i][0])
                if idsRegPull[i][1] == '*':
                        representativeId = idsRegPull[i][0]        
        # Pull out counts from htseq output
        tempCount = 0
        for trans in transIds:
                tempCount += int(htseqDict[trans])
        # Format output row
        outputText.append(representativeId + '\t' + str(tempCount))
        # Indication of progress
        ongoingCount += 1
        if ongoingCount % 1000000 == 0:
                print(str(ongoingCount) + ' clusters have been processed after ' + str(time.time() - start_time) + ' seconds')

# Output file
output = '\n'.join(outputText)
outFile = open(outName, 'w')
outFile.write(output)
outFile.close()
print('The final cluster has been processed after ' + str(time.time() - start_time) + ' seconds and output file saved.')
