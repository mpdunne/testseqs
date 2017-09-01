#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random
import sys
import string

aas=['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']

def getSeqs(errorprob, gapprob, length, numseqs, source=aas):
     ranseq = "M" + string.join([random.choice(source) for i in range(0,length)], "")
     # randomly permute the elements up to a certain probability
     seqs = [string.join([(a if random.random() < (1-errorprob) else (random.choice(aas) if random.random() < (1-gapprob) else "-")) for a in ranseq],"") for i in range(0,numseqs)]
     for i in seqs:
          print ">sequence_"+string.join([random.choice(["A","B","C","D","E","F","G","H","I"]) for p in range(8)],"")
          print i

args = sys.argv[1:]
errorprob = float(args[0])
gapprob = float(args[1])
length = int(args[2])
numseqs = int(args[3])

getSeqs(errorprob, gapprob, length, numseqs, source=aas)
