#!/usr/bin/python

##THIS SCRIPT IS MEANT TO GO THROUGH THE .LOCI FILE AND COUNT PRESENCE AND ABSENCE OF EACH LOCUS, THEN PUT INTO A TABLE

import sys
import itertools
from Bio import SeqIO
import numpy as np

infile = sys.argv[1]
statfile = sys.argv[2]
outfile = sys.argv[3]
sample_outfile = sys.argv[4]
data = open(infile).readlines()
stat = open(statfile).readlines()
samp = open(sample_outfile,'w')

samples = []
for line in stat[8:]:
    if line == '\n':
        break
    L = line.lstrip().rstrip().split("\t")
    L[0] = L[0].replace(" ","")
    samples.append(L[0])

four =  stat[4].lstrip().rstrip().split(" ")
num_loc = int(four[0])
dim = (len(samples),num_loc)        
APM = np.zeros(dim,dtype=int) ##APM is absence-presence matrix

counter = 0    
for line in data:
    a = line.lstrip().rstrip().split(" ")
    if '>' in a[0]:
        name = a[0].replace(">","")
        row = samples.index(name)
        APM[row][counter] = 1
    elif '//' in a[0]:
        counter = counter + 1

np.savetxt(outfile, APM, delimiter=' ', fmt='%.1s')
for item in samples:
    print >> samp, item

samp.close()
