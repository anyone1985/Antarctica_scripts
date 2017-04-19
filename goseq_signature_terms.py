#! python3
# goseq_signature_terms.py
# This script analyses a directory containing tab-delimited text files of all pairwise comparisons of treatments. These files contain two columns,
# Column 1: Terms
# Column 2: FDR values
# with or without a header. A user specified significance value can be used as a cut-off.
# The result are text files which are intended to identify "signature" terms unique to each site.
# These signature terms are terms that
#(1) are consistently up or downregulated when compared to all other sites and
#(2) are unique to a site i.e., which terms were not also found in (1) for any other site.
#
# Note that this script used to be two separate scripts, thus this script is set up in a step1-step2 system. Also, it means it's probably not written as efficiently or as transparently as it could be. Sorry.

# Import packages needed and define function
import os, argparse, re

def pval_associate(commonTerms, pvalDict, group, significance):
        # Re-associate P-values back to our results
        associatedList = []
        for term in commonTerms:
                # Filter P-values if user specified
                groupPval = pvalDict[term]/len(group)
                if significance != 'unspecified':
                        if float(groupPval) > float(significance):
                                continue
                associatedList.append(term + '\t' + str(groupPval))
        # Return sorted list
        sortedList = sorted(associatedList, key = lambda x: float(x.split('\t')[1]))
        return sortedList

#### USER INPUT SECTION

usage = """%(prog)s reads in all files not ending with .py extension in the current directory
 and summarises the results of all pairwise comparisons for a sample into a list of "signature" terms.
 File format must be two columns with or without a header, where column 1 contains the functional terms, and column 2 contains the P-values matching these terms. 
Note that file names must be in a format like this: SITE1.vs.SITE2.UP.txt. The spacer string ('.vs.') can be user specified and the ending string ('.UP.txt') must contain a '.' before either 'UP' or 'DOWN' to indicate the direction of comparison.
Also note that this script assumes that naming conventions are the same as how edgeR performs pairwise comparison. i.e., SITE1.vs.SITE2 is the same order of sample names used when running exactTest(data, pair=c(SITE1,SITE2) as an example.
"""

p = argparse.ArgumentParser(description=usage)
#p.add_argument("-i", "--i", "-input", "--input", dest="xml",
#                  help="k-SLAM xml file")
p.add_argument("-e", "--e", "-e", "--e", dest="significance", default = "unspecified",
                  help="E-value cut-off. If unspecified, no cut-off is used. Can be entered as an integer or in scientific notation (e.g., 1e-5)")
p.add_argument("-d", "--d", "-direction", "--direction", dest="direction", choices = ["up", "down", "both", "UP", "DOWN", "BOTH"], default = "both",
                  help="Specify if you want to only summarise up files, down files, or both.")
p.add_argument("-s", "--s", "-spacer", "--spacer", dest="spacer", default = ".vs.",
                  help="Specify the string that separates the two sample names. By default this script assumes a format like 'SITE1.vs.SITE2' in the same order as edgeR uses when performing pairwise comparison.")
args = p.parse_args()
significance = args.significance
direction = args.direction
spacer = args.spacer

#### DATA PREPARATION

# Create output directory if needed
outputDir = 'signature_terms'
if not os.path.isdir(os.getcwd() + '//' + outputDir):
        os.mkdir(os.getcwd() + '//' + outputDir)

# Sort out file naming [since the spacer string is specified, if the files are named appropriately this should save the user from manually specifying their sample names]
fileNames = []
for entry in os.listdir('.'):
        if not entry.endswith('.py') and os.path.isfile(entry):
                fileNames.append(entry)

spaceRegex = re.compile(r'(.+?)' + spacer + r'(.+?)\.')
samples = []
for name in fileNames:
        regHit = spaceRegex.search(name)
        samples.extend([regHit.group(1),regHit.group(2)])
samples = list(set(samples))            # Non-redundant sample name list

# Organise samples variable to be in proper descending order [i.e., one sample should be on the left side of the spacer N-1 times, then the next should be N-2 times, etc.]
        ## This is not necessary for the script to function, but I left it in because it ordered the samples in the same way as the edgeR output does, so it should look neater when printed to screen below
for i in range(len(samples)):
        samples[i] = [samples[i],0]       # Make each entry a list so we can append a numeric value to it
for i in range(len(samples)):
        for name in fileNames:
                if name.startswith(samples[i][0] + spacer):
                        samples[i][1] += 1
samples = list(reversed(sorted(samples, key = lambda x: x[1])))
templist = []
for sample in samples:
        templist.append(sample[0])
samples = templist

# Let user know what samples the script has detected to know if something has gone wrong
print('It looks like the samples you\'re using are: ')
for i in range(0, len(samples), 5):
        print(samples[i:i+5])
print('If this doesn\'t look right, make sure your file names are appropriately in a format identical or similar to \'SITE1.vs.SITE2.UP.txt\'')

#### CORE LOOP

### STEP 1 - Find terms that are consistently up or downregulated when compared to all other sites
step1Up = []
step1Down = []
# Loop through files in directory
for i in range(0, len(samples)):
        ### PULL OUT FILE NAMES FOR PROCESSING IN THIS SPECIFIC LOOP
        # Get all files with the current sample name in it
        currSampleFiles = []
        for file in fileNames:
                if samples[i] + '.' in file:
                        currSampleFiles.append(file) 
        # Sanity check
        if len(currSampleFiles) != (len(samples)-1)*2:
                print('There\'s either files missing in the current directory, or too many files. My calculations suggest there should be ' + str(int((len(samples)-1)*len(samples))) + ' files in this directory (excluding the python script). Check to see if there are any additional or missing files, otherwise, make sure file names are in the correct format.')
                quit()
        # Figure out which files are up and down in the sample
        upFileNames = []
        downFileNames = []
        for fname in currSampleFiles:
                if fname[0:len(samples[i])] == samples[i] and "DOWN" in fname:    # If it's (for example) SITE1 vs. SITE2 DOWN, that means the terms are DOWN in SITE2 when compared to SITE1. Or, in other words, UP in SITE1 when compared to SITE2. This is according to edgeR naming conventions.
                        upFileNames.append(fname)
                elif fname[0:len(samples[i])] == samples[i] and "UP" in fname:    # Inverse of above
                        downFileNames.append(fname)
                elif fname[0:len(samples[i])] != samples[i] and "DOWN" in fname:  # If it's (for example) SITE2 vs. SITE1 DOWN, that means the terms are DOWN in AR when compared to DW.
                        downFileNames.append(fname)
                elif fname[0:len(samples[i])] != samples[i] and "UP" in fname:    # Inverse of above
                        upFileNames.append(fname)
                else:
                        # Sanity check
                        print('There is something wrong with the file naming. Make sure all files are in a format like \'SITE1.vs.SITE2.UP.txt\'')
                        quit()
        ### PARSE FILES AND FIND UNIQUE TERMS
        # Sorting for if we're doing up, down, or both [direction value]
        if direction.lower() == 'up':
                combinedNames = [upFileNames]
        elif direction.lower() == 'down':
                combinedNames = [downFileNames]
        else:
                combinedNames = [upFileNames, downFileNames]
        pvalDict = {}                                   # This value allows us to reassociate Pvals after unique term finding
        # Loop into up and down category files
        for y in range(len(combinedNames)):
                fileContentsList = []
                pvalDictionaryList = []
                group = combinedNames[y]
                # Read and parse files into a list
                for file in group:
                        contents = open(file, 'r').read()
                        contents = contents.split('\n')
                        # Sort out potential header
                        headerCheck = contents[0].split('\t')
                        try:
                                float(headerCheck[1])   # If this succeeds, then this value is a P-value and not a header
                                header = None
                        except:
                                header = contents[0]
                                del contents[0]
                        # Remove blank lines at end, if present
                        while contents[-1] == '':
                                del contents[-1]
                        # Pull out columns from each line of the file contents
                        tempFileContents = []
                        tempDictionary = {}
                        for line in contents:
                                line = line.split('\t')
                                term = line[0]
                                pval = line[1]
                                # Drop the annotation terms into a temporary value to later add to fileContentsList for processing
                                tempFileContents.append(term)
                                # Make dictionary to re-associate P-values after sorting process [note: adds P-values here, will divide to obtain average later]
                                if term not in pvalDict:
                                        pvalDict[term] = float(pval)
                                else:
                                        pvalDict[term] = pvalDict[term] + float(pval)
                        fileContentsList.append(tempFileContents)
                        
                # Identify common terms
                for z in range(0, len(fileContentsList)):
                        if z == 0:
                                commonTerms = set(fileContentsList[z])
                        else:
                                commonTerms = commonTerms & set(fileContentsList[z])
                                if len(commonTerms) == 0:
                                        break
                # Create step1 lists
                if direction.lower() == 'up':
                        if len(commonTerms) == 0:
                                outName = samples[i] + '.UP.vs.ALL.txt'
                                step1Up.append(['None'])
                        else:
                                # Re-associate P-values back to our results
                                sortedList = pval_associate(commonTerms, pvalDict, group, significance)
                                step1Up.append(sortedList)
                elif direction.lower() == 'down':
                        if len(commonTerms) == 0:
                                outName = samples[i] + '.DOWN.vs.ALL.txt'
                                step1Down.append(['None'])
                        else:
                                # Re-associate P-values back to our results
                                sortedList = pval_associate(commonTerms, pvalDict, group, significance)
                                step1Down.append(sortedList)
                else:
                        if y == 0:
                                if len(commonTerms) == 0:
                                        outName = samples[i] + '.UP.vs.ALL.txt'
                                        step1Up.append(['None'])
                                else:
                                        # Re-associate P-values back to our results
                                        sortedList = pval_associate(commonTerms, pvalDict, group, significance)
                                        step1Up.append(sortedList)
                        else:
                                if len(commonTerms) == 0:
                                        outName = samples[i] + '.DOWN.vs.ALL.txt'
                                        step1Down.append(['None'])
                                else:
                                        # Re-associate P-values back to our results
                                        sortedList = pval_associate(commonTerms, pvalDict, group, significance)
                                        step1Down.append(sortedList)

# Pre-step 2 setup
step2List = []
if direction.lower() == 'up':
        step2List = [step1Up]
        wordlist = ['UP']
elif direction.lower() == 'down':
        step2List = [step1Down]
        wordlist = ['DOWN']
else:
        step2List = [step1Up, step1Down]
        wordlist = ['UP', 'DOWN']

print('')

#### STEP 2 - Find terms that are unique to a site i.e., terms that were not also found in (1) for any other site.
for i in range(0, len(samples)):
        for n in range(len(wordlist)):
                # Get current list
                currFile = step2List[n][i]
                # Handle 'None' files
                if currFile == ['None']:
                        outName = samples[i] + '.' + wordlist[n] + '.None.txt'
                        print('Found nothing for this group (' + samples[i] + ' ' + wordlist[n] + ')')
                        # Create quick output file
                        outputFile = open(os.getcwd() + '//' + outputDir + '//' + outName, 'w')
                        outputFile.write('None')
                        outputFile.close()
                        continue
                # Map P-values for later reassociation
                currPvalDict = {}
                temp = []
                for line in currFile:
                        line = line.split('\t')
                        currPvalDict[line[0]] = line[1]
                        temp.append(line[0])                    # Do this to hijack script written for slightly different purpose - use this to remove the P-val and just have the GO term
                currFile = temp
                # Load in remaining files
                fileContentsList = []
                for x in range(0, len(samples)):
                        if x == i:
                                continue
                        else:
                                contents = step2List[n][x]
                                temp = []
                                for line in contents:
                                        temp.append(line.split('\t')[0])                        # As above, doing this to hijack old script
                                contents = temp
                                fileContentsList.append(contents)
                # Identify unique terms
                for y in range(0, len(fileContentsList)):
                        currContents = fileContentsList[y]
                        for entry in currContents:
                                if entry in currFile:
                                        del currFile[currFile.index(entry)]
                # Check remaining terms
                if len(currFile) == 0:
                        # Sort out output name
                        outName = samples[i] + '.' + wordlist[n] + '.None.txt'
                        print('Found nothing for this group (' + samples[i] + ' ' + wordlist[n] + ')')
                        # Create quick output file
                        outputFile = open(os.getcwd() + '//' + outputDir + '//' + outName, 'w')
                        outputFile.write('None')
                        outputFile.close()
                # Format output if we do have common terms
                else:
                        # Reassociate P-values
                        associatedTerms = []
                        for line in currFile:
                                associatedTerms.append(line.split('\t')[0] + '\t' + currPvalDict[line.split('\t')[0]])
                        # Create file now
                        outName = samples[i] + '.' + wordlist[n] + '.vs.ALL.txt'
                        outputText = '\n'.join(associatedTerms)
                        outputFile = open(os.getcwd() + '//' + outputDir + '//' + outName, 'w')
                        outputFile.write(outputText)
                        outputFile.close() 

