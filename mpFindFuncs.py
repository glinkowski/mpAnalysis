# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# functions used to Find Paths in the processed network
#
# These functions were created to aid in the reading of a
#	sample, creation of random samples, and calculation
#	of statistics pertaining to the metapaths for a sample.
#   
#
# Functions provided:
#	readFileAsList(fname)
#	readFileAsIndexDict(fname)
#	readSampleFiles(fname, up, down)
#	readKeyFile(path, name)
#	readGenesFile(path, name)
#	checkGenesInNetwork(path, name, geneList)
#	convertToIndices(names, iDict)
#	getPathMatrix(mpTuple, path, name)
#	getPathCountOne(sample, matrix)
#	getPathCountList(samples, matrix)
#	getPathMeans(rSamples, matrix)
#	calculateStatistics(sample, rSamples, mpDict,
#			path, name)
#	nameOutputFile(path, name)
#	writeOutputOneSample(path, name, nName, sName,
#			mpDict, counts, means, stdevs, scores, leftout)
# ---------------------------------------------------------

import os.path
from os import listdir
import sys
import numpy as np
import random



####### ####### ####### ####### 
# PARAMETERS

# Data type used when loading edge file to memory:
#nodeDT = np.dtype('a30')
# Length to pad the matrix file names:
keyZPad = 5
# Length to pad the output file names:
fnZPad = 3
# Data delimiter to use in the output file:
textDelim = '\t'
# Whether to save uncompressed text version of matrix:
#saveText = True		# (useful for error-checking)

####### ####### ####### ####### 



######## ######## ######## ######## 
# Function: Read in a file as a line-by-line list of items
# Input ----
#   fname, str: path + name of the the sample files
# Returns ----
#	fItems, str list: ordered list of items from file
def readFileAsList(fname) :

	# ERROR CHECK: verify file exists
	if not os.path.isfile(fname) :
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	# The list of items to return
	fItems = list()

	# Read in from the file
	fn = open(fname, "rb")
	for line in fn :
		fItems.append( line.rstrip() )
	#end loop
	fn.close()

	fItems.sort()
	return fItems
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Read in the gene file. File is an ordered list
#	of genes, where the row number (starting at zero)
#	corresonds to the index in the matrix/list/etc where
#	the gene can be found.
# Input ----
#	fname, str: path & name to keep file
# Returns ----
#	iDict, dict: 
#		key, str: gene name as read from file
#		value, int: index to corresponding array
def readFileAsIndexDict(fname) :

	# ERROR CHECK: verify file exists
	if not os.path.isfile(fname) :
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	# Build the dictionary from the text file
	iDict = dict()
	gf = open(fname, "rb")
	index = 0
	for line in gf :
		gene = line.rstrip()	# remove "\n"
		iDict[gene] = int(index)
		index += 1
	#end loop

	return iDict
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Read in the dataset from a samplename
#	Check for variants: ".txt", "_UP.txt", "_DN.txt"
# Input ----
#   fname, str: path + name of the the sample files
#   up, down: boolean -- only read the _UP or _DN file if true
# Returns ----
#	sNodes, str list: ordered list of names from file(s)
def readSampleFiles(fname, up, down) :

	# The list of items to return
	sNodes = list()

	# Flag indicates a file existed and was read
	exists = False

	# First look for the file as named (no _UP or _DN)
	if os.path.isfile(fname + ".txt") :
		temp = readFileAsList(fname + ".txt")
		sNodes.extend(temp)
		exists = True
	#end if

	# Look for the _DN file
	if down and os.path.isfile(fname + "_DN.txt") :
		temp = readFileAsList(fname + "_DN.txt")
		sNodes.extend(temp)
		exists = True
	#end if

	if up and os.path.isfile(fname + "_UP.txt") :
		temp = readFileAsList(fname + "_UP.txt")
		sNodes.extend(temp)
		exists = True
	#end if

	# Alert user if nothing was read in
	if not exists :
		print "WARNING: no file found: {}".format(fname)

	# Do NOT return duplicates
	uNodes = np.unique(sNodes) # sorted list of unique items
	return uNodes
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Read in the key.txt file regarding the 
#	metapath matrices
# Input ----
#	path, str: path to the network files
#	name, str: name of the network to use
# Returns ----
#	keyDict, dict
#		key, str: name of metapath
#		value, tuple: int is matrix/file ID number
#			bool where True means use matrix transpose
def readKeyFile(path, name) :

	fname = path + name + "_MetaPaths/key.txt"

	# ERROR CHECK: verify file exists
	if not os.path.isfile(fname) :
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	# The item to return
	keyDict = dict()

	# Read in the file
	fk = open(fname, "rb")
	firstline = True
	for line in fk :

		# skip the first line
		if firstline :
			firstline = False
			continue
		#end if

		# separate the values
		line = line.rstrip()
		lk = line.split('\t')
		lv = lk[0].split(',')

		transpose = False
		if lv[1] == "t" :
			transpose = True
		#end if

		# add to the dict
		keyDict[lk[1]] = [int(lv[0]), transpose]
	#end loop
	fk.close()

	return keyDict
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Read in the genes.txt file containing the 
#	gene-name headers to the metapath matrices
# Input ----
#	path, str: path to the network files
#	name, str: name of the network to use
# Returns ----
#	gDict, dict
#		key, str: name of gene
#		value, int: row/col index for that gene
def readGenesFile(path, name) :

	fname = path + name + "_MetaPaths/genes.txt"

	# The item to return
	gDict = readFileAsIndexDict(fname)

	return gDict
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Take the list of genes from the sample, return
#	a list of genes in the network, and left-out genes
# Input ----
#	path, str: path to the network files
#	name, str: name of the network to use
#	geneList, str list: genes in the sample
# Returns ----
#	inList & outList, str list: which genes from the
#		given list are found in network & which aren't
def checkGenesInNetwork(path, name, geneList) :

	fname = path + name + "_MetaPaths/genes.txt"

	# ERROR CHECK: verify file exists
	if not os.path.isfile(fname) :
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	# Read in the genes from the file
	geneSet = set()
	fg = open(fname, "rb")
	for line in fg :
		line = line.rstrip()
		geneSet.add(line)
	#end loop
	fg.close()

	# The items to return
	inList = list()
	outList = list()

	# Sift through the sample
	for gene in geneList :
		if gene in geneSet :
			inList.append(gene)
		else :
			outList.append(gene)
		#end if
	#end loop

	inList.sort()
	outList.sort()
	return inList, outList
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Take the list of genes from the sample, return
#	a list of genes in the network, and left-out genes
# Input ----
#	numRows, int: number of random samples to create
#	numCols, int: number of items in each random sample
#	numItems, int: number of items in the source
#		(ie: the list from which the rand sample is pulled)
# Returns ----
#	rSamples, int array: 
def createRandomSamplesArray(numRows, numCols, numItems) :

	# The item to return
	rSamples = np.empty([numRows, numCols], dtype=np.int32)

	# Choose the random samples
	for r in range(0, numRows) :
		rSamples[r,:] = random.sample(xrange(numItems),
			numCols)
	#end loop

	return rSamples
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Given a list of names and appropriate dict,
#	convert to corresponding index values
# Input ----
#	names, str list: names of items to convert
#	iDict, dict
#		key, str: names of items
#		value, int: corresponding index values
# Returns ----
#	indices, int list: index values of the input items
def convertToIndices(names, iDict) :

	# The item to return
	indices = list()

	for name in names :
		indices.append(iDict[name])
	#end loop

	return indices
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Find the number of paths of this type joining
#	the nodes in the sample
# Input ----
#	mpTuple [int, bool]: indicates which matrix file to use
#	path, str: path to the network files
#	name, str: name of the network to use
# Returns ----
#	matrix, int array: num paths between node pairs
def getPathMatrix(mpTuple, path, name) :

	zpad = keyZPad
	fname = (path + name + "_MetaPaths/" +
		"{}.npy".format(str(mpTuple[0]).zfill(zpad)) )

	# ERROR CHECK: verify file exists
	if not os.path.isfile(fname) :
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	# Load the matrix
	matrix = np.load(fname)

	# Convert to transpose if flag==True
	if mpTuple[1] :
		return np.transpose(matrix)
	else :
		return matrix
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Find the number of paths of this type joining
#	the nodes in the sample
# Input ----
#	sample, int list: indices of the nodes in the sample
#	matrix, int array: num paths between node pairs
# Returns ----
#	count, int: num paths of this type in this sample
def getPathCountOne(sample, matrix) :

#TODO: If there is no difference in counts between a path
#	and it's inverse, then the inverse path need not be
#	calculated.
#NOTE: If a particular path is chosen for prediction, then
#	its inverse should be evaluated as well. (simulaneously)

	# Count the paths in the sample
	count = 0
	for a in sample :
		for b in sample :
			if a != b :		# skip if a==b
				count += matrix[a,b]
	#end loop

	return count
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Find the number of paths of this type joining
#	the nodes in the sample
# Input ----
#	samples, int 2D array: indices of the nodes in the sample
#		rows: individual samples
#	matrix, int array: num paths between node pairs
# Returns ----
#	counts, int list: num paths of this type in each sample
def getPathCountList(samples, matrix) :

	counts = list()
	for i in range(0, samples.shape[0]) :
		# Count the paths in the sample
		count = 0
		for a in samples[i,:] :
			for b in samples[i,:] :
				if a != b :		# skip if a==b
					count += matrix[a,b]
		#end loop
		counts.append(count)
	#end loop

	return counts
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: For a set of samples, find the mean and standard
#	devation of path counts for a specific path
# Input ----
#	samples, int 2D array: indices of the nodes in the sample
#		rows: individual samples
#	matrix, int array: num paths between node pairs
# Returns ----
#	mean, float list: mean count of paths for each sample
#	stdev, float list: standard deviation of counts per samp
def getPathMeans(rSamples, matrix) :

	# Get the path counts for each random sample
	counts = getPathCountList(rSamples, matrix)

	# Return the mean  and standard deviation
	mean = np.mean(counts)
	stdev = np.std(counts)
	return mean, stdev
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Calculate the statistics necessary to find
#	outstanding paths within the sample
# Input ----
#	samples, int 2D array: indices of the nodes in the sample
#		rows: individual samples
#	mpDict, {str: [int, bool]} dict:
#		key, str - name of the metapath
#		value, [int, bool] - which matrix file to use, and 
#			whether to use the transpose (inverse path)
#	path, str: path to the network files
#	name, str: name of the network to use
# Returns ----
#	sCount, int list: num paths of this type in this sample
#	rMeans, float list: mean count of path in rand samples
#	rStDev, float list: standard deviation of rand samples
#	zScore, float list: z-Score for each path type
def calculateStatistics(sample, rSamples, mpDict,
	path, name) :

	# The items to return
	sCount = list()		# num of each path in given sample
	rMeans = list()		# mean count of e. p. in rand samples
	rStDev = list()		# stand dev of e. p. in rand samples
	zScore = list()		# z-Score of e. p. in rand samples

	# An iterable list of metapaths
	mpList = mpDict.keys()
	mpList.sort()

	for mp in mpList :

		matrix = getPathMatrix(mpDict[mp], path, name)

		tCount = getPathCountOne(sample, matrix)
		sCount.append(tCount)

		tMeans, tStDev = getPathMeans(rSamples, matrix)
		rMeans.append(tMeans)
		rStDev.append(tStDev)

		tScore = (tCount - tMeans) / tStDev
		zScore.append(tScore)
	#end loop

	return sCount, rMeans, rStDev, zScore
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: choose an unused name for the output file
# Input ----
#	path, str: path to the network files
#	name, str: name of the network to use
# Returns ----
#	fname, str: name of output file (without path)
def nameOutputFile(path, name) :

	# ERROR CHECK: verify directory exists
	if not os.path.isdir(path) :
		print ( "ERROR: Specified path doesn't exist:" +
			" {}".format(path) )
		sys.exit()
	#end if

	zpad = fnZPad

	# Set of all files in the directory
	fileSet = set(os.listdir(path))

	# increment file name until an unused one is found
	num = int(0)
	fname = name + "-{}.txt".format(str(num).zfill(zpad))
	while fname in fileSet :
		num += 1
		fname = name + "-{}.txt".format(str(num).zfill(zpad))
	#end loop

	return fname
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: For a single sample, write the data to a file
# Input ----
#	path, str: path to the output files
#	name, str: name to use in creating output file
#	nName, str: name of the network used
#	sName, str: name of the sample used
#	mpDict, {str: [int, bool]} dict:
#		key, str - name of the metapath
#		value, [int, bool] - which matrix file to use, and 
#			whether to use the transpose (inverse path)
#	counts, int list: num paths of this type in this sample
#	means, float list: mean count of path in rand samples
#	stdevs, float list: standard deviation of rand samples
#	scores, float list: z-Score for each path type
#	leftout, str list: list of genes in original sample
#		which don't appear in the network
# Returns ----
#	fname, str: name of output file (without path)
def writeOutputOneSample(path, name, nName, sName,
	mpDict, counts, means, stdevs, scores, leftout) :

	delim = textDelim

	# Name the output file
	fname = nameOutputFile(path, name)

	# Get the ordered list of metapaths
	mpList = mpDict.keys()
	mpList.sort()

	# Write the file header
	fn = open(path+fname, "wb")
	fn.write("The occurrence of metapaths with in a sample ...\n")
	fn.write("network:{}{}\n".format(delim, nName))
	fn.write("sample:{}{}\n".format(delim, sName))
	fn.write("\n")

	# Write the body of the file
	fn.write("sample{}random samples\n".format(delim))
	fn.write("count{0}mean{0}std dev{0}".format(delim) +
		"z-score{0}path name\n".format(delim))
	i = 0
	for mp in mpList :
		fn.write("{1}{0}{2}{0}{3}{0}{4}{0}{5}\n".format(delim,
			counts[i], means[i], stdevs[i], scores[i], mp) )
		i += 1
	#end loop

	# Write the file footer
	fn.write("\n")
	fn.write("Sample genes not found in the network:" +
		" {}".format(len(leftout)) )
	for gene in leftout :
		fn.write("\n{}".format(gene) )
	#end loop
	fn.close()

	return fname
#end def ######## ######## ######## 