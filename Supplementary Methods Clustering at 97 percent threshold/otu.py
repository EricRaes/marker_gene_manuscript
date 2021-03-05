#!/usr/bin/env python

"""
Make a list of seqs from a zotu table, with size annotation.
Takes the input file as an argument and uses that file to create a tab delimited file comprising the 
sequence and the sequences plus the size annotation. Output is to stdout
"""

from __future__ import print_function

import sys
import numpy as np

bytes = 0

for i,line in enumerate(open(sys.argv[1])):
    bytes += len(line)
    if i == 0: continue
    p = line.index('\t')
    otu = line[0:p]
    a = np.fromstring(line[p+1:], dtype=float, sep='\t')
    print(otu, otu + ';size=' + str(np.sum(a)), sep='\t')
    if i%10000 == 0: print(i, bytes, file=sys.stderr)
