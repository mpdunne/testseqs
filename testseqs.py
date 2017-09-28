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

def cry():
	print "Invalid input: Please use format\n\n\t\t\ttestseqs p_err p_gap len numseqs\n\n\t\twhere\n\
									\n\t\t\tp_err   -- is the error probability \
									\n\t\t\tp_gap   -- is the probability that a given error is a gap \
									\n\t\t\tlength  -- is the length of the alignment\
									\n\t\t\tnumseqs -- is the number of sequences\n"
	sys.exit()

args = sys.argv[1:]
try:
	errorprob = float(args[0])
	gapprob = float(args[1])
	length = int(args[2])
	numseqs = int(args[3])
except:
	cry()

if not (0 <= errorprob <= 1 and 0 <= gapprob <= 1 and length > 0 and numseqs > 0): cry()

getSeqs(errorprob, gapprob, length, numseqs, source=aas)
