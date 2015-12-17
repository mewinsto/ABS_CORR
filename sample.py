#!/usr/bin/python

import sys
import itertools
import numpy as np

infile = sys.argv[1]
outfile = sys.argv[2]
stat = open(infile).readlines()
samp_out = open(outfile, 'w')

samples = []
for line in stat[8:]:
    if line == '\n':
        break
    L = line.lstrip().rstrip().split("\t")
    L[0] = L[0].replace(" ","")
    samples.append(L[0])

for i in range(len(samples)):
    samp_out.write(samples[i])
    samp_out.write("\n")

samp_out.close()
