#! python3
# blast_hit_counter.py
# Simple script to count the number of sequences that obtain a
# hit more significant than a specified cut-off

import argparse
from itertools import groupby

usage = """%(prog)s reads BLAST-tab files and counts how many sequences
have a hit greater than a certain E-value
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="file",
                  help="BLAST-tab file")
p.add_argument("-e", "--e", "--evalue", "-evalue", type=float, dest="evalue",
                   help="E-value cut-off", default=1e-05)
args = p.parse_args()

hitCount = 0
grouper = lambda x: x.split('\t')[0]
with open(args.file) as fileIn:
        for key, value in groupby(fileIn, grouper):
                bestHit = []
                for entry in value:     # We do this process in case the file isn't sorted; this is the case for MMseqs2
                        line = entry.rstrip('\n').rstrip('\r').split('\t')
                        if bestHit == []:
                                bestHit = line
                        else:
                                # Check E-value
                                if float(line[10]) < float(bestHit[10]):
                                        bestHit = line
                                # If E-value is equivalent, check bit score
                                elif float(line[10]) == float(bestHit[10]):
                                        if float(line[11]) > float(bestHit[11]):
                                                bestHit = line
                                        # If bit score is equivalent, check alignment length [this is arbitrary, if E-value and bit score are identical, then the possible scenarios are that one hit is longer with lower identity, and the other hit is shorter with higher identity. I'm inclined to think that the first hit is better than the second if E-value/bit score are equivalent...]
                                        elif float(line[11]) == float(bestHit[11]):
                                                if int(line[3]) > int(bestHit[3]):
                                                        bestHit = line
                # Process line to format it for output
                if float(bestHit[10]) > args.evalue:
                        continue
                else:
                        hitCount += 1

# Print output
print(args.file + ' num hits more significant than ' + str(args.evalue) + ': ' + str(hitCount))
