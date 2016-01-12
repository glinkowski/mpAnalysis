# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# functions used to Find Paths in the processed network
#
# These functions were created to aid in the
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
# ---------------------------------------------------------

import os.path
import sys
import numpy as np
import random



####### ####### ####### ####### 
# PARAMETERS

# Data type used when loading edge file to memory:
#nodeDT = np.dtype('a30')
# Length to pad the matrix file names:
keyZPad = 5
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

#	# maxVal is maximum index value to return
#	maxVal = numItems - 1

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
#	samples, int 2D array: indices of the nodes in the sample
#		rows: individual samples
#	mpTuple [int, bool]: indicates which matrix file to use
#	path, str: path to the network files
#	name, str: name of the network to use
# Returns ----
#	count, int: num paths of this type in this sample
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
#	mpTuple [int, bool]: indicates which matrix file to use
#	path, str: path to the network files
#	name, str: name of the network to use
# Returns ----
#	count, int: num paths of this type in this sample
def getPathCountOne(sample, matrix) :

#	zpad = keyZPad
#	fname = (path + name + "_MetaPaths/" +
#		"{}.npy".format(str(mpTuple[0]).zfill(zpad)) )
#
#	# ERROR CHECK: verify file exists
#	if not os.path.isfile(fname) :
#		print ( "ERROR: Specified file doesn't exist:" +
#			" {}".format(fname) )
#		sys.exit()
#	#end if
#
#	# Load the matrix
#	matrix = np.load(fname)
##	print matrix

#	# Load the matrix
#	matrix = getPathMatrix(mpTuple, path, name)
#
#	# Count the paths in the sample
#	count = 0
#	# If flag==True, use the matrix transpose
##TODO: Verify this: the path counts should be the same.
#	if mpTuple[1] == False :
#		for a in sample :
#			for b in sample :
#				if a != b :		# skip if a==b
#					count += matrix[a,b]
#	else :
#		for a in sample :
#			for b in sample :
#				if a != b :		# skip if a==b
#					count += matrix[b,a]
#	#end loop

#TODO: If there is no difference in counts between a path
#	and it's inverse, then the inverse path need not be
#	calculated.
#NOTE: If a particular path is chosen for prediction, then
#	its inverse should be evaluated as well.


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
#	mpTuple [int, bool]: indicates which matrix file to use
#	path, str: path to the network files
#	name, str: name of the network to use
# Returns ----
#	count, int: num paths of this type in this sample
def getPathCountList(samples, matrix) :

#	zpad = keyZPad
#	fname = (path + name + "_MetaPaths/" +
#		"{}.npy".format(str(mpTuple[0]).zfill(zpad)) )
#
#	# ERROR CHECK: verify file exists
#	if not os.path.isfile(fname) :
#		print ( "ERROR: Specified file doesn't exist:" +
#			" {}".format(fname) )
#		sys.exit()
#	#end if
#
#	# Load the matrix
#	matrix = np.load(fname)
##	print matrix

#	# Load the matrix
#	matrix = getPathMatrix(mpTuple, path, name)

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
#	mpTuple [int, bool]: indicates which matrix file to use
#	path, str: path to the network files
#	name, str: name of the network to use
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
#	mpTuple [int, bool]: indicates which matrix file to use
#	path, str: path to the network files
#	name, str: name of the network to use
# Returns ----
#	count, int: num paths of this type in this sample
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
