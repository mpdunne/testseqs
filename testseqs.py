#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2017 Michael Dunne
#
# Testseqs
#
# For any enquiries send an email to Michael Dunne
# michael.dunne@worc.ox.ac.uk

####################################
###### Safely import packages ######
####################################

errors   = []
libnames = ["random", "sys", "string", "argparse", "math", "csv"]

for libname in libnames:
    try:
	lib = __import__(libname)
    except ImportError as e:
	errors.append(e)
    else:
	globals()[libname] = lib

try:
	from bisect import bisect
except ImportError as e:
	errors.append(e)

try:
	from ete3 import Tree
except ImportError as e:
	errors.append(e)


try:
	from collections import Counter
except ImportError as e:
	errors.append(e)

try:
	from Bio import SeqIO
except ImportError as e:
	errors.append(e)

if errors:
	print("Missing modules :(\nThe following module errors need to be resolved before running Testseqs:")
	for error in errors: print("-- " + str(error))
	sys.exit()

###################################
###### Function definitions  ######
###################################

aas=['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']

def readCsv(path_file, ignoreBlank = True, ignoreHash = False):
	with open(path_file, "r") as f:
		data = [line for line in csv.reader(f, delimiter="\t") if \
			(not ignoreBlank or ''.join(line).strip()) and (not ignoreHash or not line[0].startswith('#'))]
	return data

def wchoice(choices):
	values, weights = zip(*choices)
	total = 0; cum_weights = []
	for w in weights:
		total += w
		cum_weights.append(total)
	x = random.random() * total
	i = bisect(cum_weights, x)
	return values[i]

def cry():
	sys.exit("Invalid input :(")

def getSeqs(ranseq, errorprob, gapprob, numseqs, weights, includestop, startaa, source=aas):
	# randomly permute the elements up to a certain probability
	seqs = [string.join([(a if random.random() < (1-errorprob) else (wchoice(weights) if random.random() < (1-gapprob) else "-")) for a in ranseq],"") for i in range(0,numseqs)]
	q = max(5, math.log(len(seqs), 10) + 1)
	for i,s in enumerate(seqs):
		if includestop: s = s[:-1] + "*"
		if startaa: s = startaa + s[1:]
		print ">sequence_" + str(i).zfill(q)
		print s
	sys.exit()

def getFWeights(p_fa, aas):
	try:
		res = ""
		for s in SeqIO.parse(p_fa, "fasta"):
			res += s
		c = Counter(res)
		aares = [(aa, c[aa]) for aa in aas]
		assert(list(set(c.keys())) == list(set(aas)))
		return aares
	except:
		sys.exit("Error: could not read FASTA file " + p_fa)

def getWeights(p_w, aas):
	try:
		w = readCsv(p_w)
		res = [(a[0], a[1]) for a in w]
		assert(list(set(res.keys())) == list(set(aas)))
		
	except:
		sys.exit("Error: Invalid or missing weights file")

def parseBool(sbool, sname):
	if sbool.lower() == "true":
		return True
	elif sbool.lower() == "false":
		return False
	else:
		sys.exit("Error: " + sname + " argument must be either true or false")

def getAas(l):
	try:
		a=list(set(l.replace(" ", "").replace(",", "")))
		assert(len(a) != 0)
		for i in a: assert(i.isalpha())
		return a
	except:
		sys.exit("Error: invalid AA values")

def getTree(p_t):
	try:
		if not p_t: return None
		return Tree(p_t, format=1)
	except:
		sys.exit("Error: Tree file " + p_t + " not found, or not in Newick format.")

####################################
########### Entry code #############
####################################

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Run TestSeqs")
	parser.add_argument("-g", "--gapprob", metavar="gapprob", dest="GP")
	parser.add_argument("-g", "--gapprob", metavar="gapprob", dest="GP")
	parser.add_argument("-e", "--errorprob", metavar="errorprob", dest="EP")
	parser.add_argument("-l", "--length", metavar="length", dest="LN")
	parser.add_argument("-n", "--numseqs", metavar="numseqs", dest="NS")
	parser.add_argument("-w", "--weights", metavar="weights", dest="WT")
	parser.add_argument("-a", "--aas", metavar="weights", dest="AA")
	#parser.add_argument("-d", "--dna", action="store_true", dest="DNA")
	#parser.add_argument("-p", "--protein", action="store_true", dest="PROT")
	parser.add_argument("-i", "--initial", metavar="initial", dest="IN")
	#parser.add_argument("-t", "--tree", metavar="tree", dest="TREE")
	parser.add_argument("-f", "--fasta", metavar="fasta", dest="FA")
	#parser.add_argument("-z", "--lossprob", metavar="loss", dest="GL")
	#parser.add_argument("-y", "--dupprob", metavar="dup", dest="GD")
	parser.add_argument("-x", "--includestop", metavar="includestop", dest="IS")
	parser.add_argument("-m", "--metstart", metavar="metstart", dest="MS")
	parser.add_argument("-c", "--altstart", metavar="altstart", dest="AS")
	#parser.add_argument("-s", "--submethod", metavar="submethod", dest="SM")
	args = parser.parse_args()
	# Process initial aa list
	if args.AA != None: aas = getAas(args.AA)
	# Get weights for initial sequence generation
	if args.WT == None and args.FA == None: weights = [(aa,1) for aa in aas]
	if args.WT != None and args.FA == None: sys.exit("Error: -w and -f arguments cannot be used together")
	if args.WT != None: weights = getWeights(args.WT, aas)
	if args.FA != None: weights = getFWeights(args.FA, aas)
	# Process initial mode
	if args.IN != None:
		initial = args.IN
		if args.LN != None: sys.exit("Error: -i and -l arguments cannot be used together")
	else:
		try:
			length = int(args.LN) if args.LN else 40
		except:
			sys.exit("Error: length argument must be an integer")
		initial = string.join([wchoice(weights) for i in range(0,length)], "")
	# Process start and stop things...
	includestop = parseBool(args.IS, "-x") if args.IS else True
	metstart    = parseBool(args.MS, "-m") if args.MS else True
	if args.AS and args.MS:
		sys.exit("Error: can't use both -m and -c options")
	if args.AS:
		startaa = str(args.AS)
		if not startaa.isalpha() and len(startaa=1): sys.exit("Error: invalid start amino acid")
	elif metstart:
		startaa = "M"
	else: startaa = ""
	minlen      = 1 + metstart + includestop + (startaa != "")
	if minlen > length: sys.exit("Error: minimum length with those options is " + str(minlen))
	# Process 
	try:
		errorprob = float(args.EP) if args.EP else 0.1
		gapprob   = float(args.GP) if args.GP else 0.2
		numseqs   = int(args.NS) if args.NS else 5
	except:
		cry()
	if not (0 <= errorprob <= 1 and 0 <= gapprob <= 1 and length > 0 and numseqs > 0): cry()
	getSeqs(initial, errorprob, gapprob, numseqs, weights, includestop, startaa, source=aas)
