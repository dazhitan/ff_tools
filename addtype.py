#!/usr/bin/env python
import sys
infile = sys.argv[1]
outfile = sys.argv[1] + 'new'
outfh = open(outfile, 'w')

oldtype, newtype = sys.argv[2:4]

newlines = []
with open(infile, 'r') as infh:
    for line in infh:
        line = line.rstrip()
        if ' %s"' % oldtype in line:
            print >> outfh, line
            newline = line.replace(' %s"' % oldtype, ' %s"' % newtype)
            newlines.append(newline)
            #print >> outfh, newline
        elif '"%s"' % oldtype in line:
            print >> outfh, line
            newline = line.replace('"%s"' % oldtype, '"%s"' % newtype)
            newlines.append(newline)
            #print >> outfh, newline
        elif '"%s ' % oldtype in line:
            print >> outfh, line
            newline = line.replace('"%s ' % oldtype, '"%s ' % newtype)
            newlines.append(newline)
            #print >> outfh, newline
        elif ' %s ' % oldtype in line:
            print >> outfh, line
            newline = line.replace(' %s ' % oldtype, ' %s ' % newtype)
            newlines.append(newline)
            #print >> outfh, newline
        else:
            print >> outfh, line

for nl in newlines:
    print >> outfh, nl
