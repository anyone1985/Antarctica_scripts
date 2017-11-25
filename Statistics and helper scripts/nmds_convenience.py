#! python3
# nmds_convenience.py
# This script was made quickly just to provide assistance with parsing a NMDS file. 
# Use it by opening the NMDS file within Microsoft Excel or equivalent program,
# and ctrl + A to highlight all content and ctrl + C to copy this. Run the program,
# press ENTER, and sorted columns of values will be placed into your clipboard to paste
# into an environment (like Microsoft Excel).

import os, pyperclip

##### USER INPUT SECTION

print('go')
buttonPress = input()
table=pyperclip.paste()
table=table.rstrip('\r\n').split('\r\n')
header = table[0].split('\t')
del table[0]
# Pull out columns
for i in range(len(table)):
        table[i] = table[i].split('\t')

siteNames = []
# Reformat
listHold = []
for i in range(len(table)):
        sampleName = table[i][0]
        siteNames.append(sampleName)
        sampleRow = '\n'.join(table[i][1:])
        preSort = []
        for x in range(1, len(table[i])):
                preSort.append([header[x], table[i][x]])
        sortedList = list(reversed(sorted(preSort, key = lambda x: float(x[1]))))
        for y in range(len(sortedList)):
                sortedList[y] = '\t'.join(sortedList[y]) + '\t'
        listHold.append(sortedList)

out = list(zip(*listHold))
for i in range(len(out)):
        out[i] = '\t'.join(out[i])
siteNames = '\t\t\t'.join(siteNames) + '\r\n'
out = '\r\n'.join(out)
pyperclip.copy(siteNames + out)
