# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# functions used by the Metapath Analysis scripts
#
#
# Functions provided:
#	verifyDirectory(path, create)
#	verifyFile(path, name, quiet)
#	saveListToText(path, name, theList)
#	saveMatrixText(matrix, mname, mpath, integer)
#	checkGenesInNetwork(path, name, geneList)
#	partitionSample(nPath, nName, oPath, sample, percent)
#	readKeyFile(path, name)
#	readFileAsIndexDict(fname)
#	readGenesFile(path, name)
#	createRandomSamplesArray(numRows, numCols, numItems)
#	readFileAsList(fname)
#	readSampleFiles(sfile, up, down)
#	convertToIndices(names, iDict)
#	getPathMatrix(mpTuple, path, name, sizeOf)
#	getSimMatrix(mpTuple, path, name, sizeOf)
#	removeInvertedPaths(mpDict)
#	getPathCountOne(sample, matrix)
#	getPathCountList(samples, matrix)
#	getPathMeans(rSamples, matrix)
#	getPercentile(oCount, rSamples, matrix)
#	getPercentDifference(oCount, rSamples, matrix)
#	calculateStatistics(sample, rSamples, mpDict, path, name)
#	nameOutputPath(path, dirPre)
#	nameOutputFile(path, name)
#	writeOutputOneSample(path, name, nName, sName, mpDict, 
#		counts, means, stdevs, scores, percents, leftout)
#	chooseTopKPathsSimple(k, ranker, mpDict)
#	applyGroupPathSim(path, name, mpTuple, sample)
#	writeItemRanks(path, statArray, itemDict, itemIndex, colHeader)
#	createRandomSamplesBinned(path, name, sample,
#		indexDict, cutoffPercent, eType, numRandSamples)
#	writeRankedGenes(path, statArray, itemDict, itemIndex,
#		hiddenSet, cutoffs)
#	writeRankedGenes02(path, suffix, statArray, itemDict,
#		itemIndex, hiddenSet, cutoffs)
#	writeRankedPaths(path, ranker, mpDict)
#	readRankedPaths(path, name)
#	writeGenericLists(path, fname, columnList)
#	getSubDirectoryList(root)
#	getMatrixDimensions(path, name)
#	
#	
#	
# ---------------------------------------------------------

import os.path
from os import listdir
import sys
import numpy as np
import random
import gzip



######## ######## ######## ######## 
# PARAMETERS

# Seed the random functions, if desired
random.seed(42)

# File names to use when partitioning gene set:
fnKnownGenes = 'known.txt'
fnConcealedGenes = 'concealed.txt'
fnLeftOutGenes = 'ignored.txt'

## Data type used when loading edge file to memory:
#nodeDT = np.dtype('a30')
nodeDT = object

# Length to pad the matrix file names:
fnMatrixZPad = 5
# Length to pad the output file names:
fnOutputZPad = 3
# Data delimiter to use in the output file:
textDelim = '\t'
# Whether to print non-error messages within these funcs
verbose = True
# Whether to use the data-type for the matrices:
speedVsMemory = False	# True favors speed, disables dtype
# Data-type for the path matrices:
matrixDT = np.float32	#TODO: any considerations here?

######## ######## ######## ######## 




######## ######## ######## ########
# Functions to set the global library parameters
# Input ----
#	NOTE: an input value of...
#		-1 keeps the parameter unchanged
#		-2 resets parameter to default
def setParamVerbose(newVal) :
	global verbose

	if newVal :
		verbose = True
	else :
		verbose = False
	#end if

	return
#end def ######## ######## ######## 
def setParamTextDelim(newVal) :
	global textDelim

	if str(newVal) == '-1' :
		textDelim = textDelim
	elif str(newVal) == '-2' :
		textDelim = '\t'
	else :
		textDelim = str(newVal)
	#end if

	return
#end def ######## ######## ######## 
def setParamFileZeroPad(newMval, newOval) :
	global fnMatrixZPad
	global fnOutputZPad

	if newMval == -1 :
		fnMatrixZPad = fnMatrixZPad
	elif newMval == -2 :
		fnMatrixZPad = 5
	else :
		fnMatrixZPad = newMval
	#end if
	if newOval == -1 :
		fnOutputZPad = fnOutputZPad
	elif newOval == -2 :
		fnOutputZPad = 3
	else :
		fnOutputZPad = newMval
	#end if

	return
#end def ######## ######## ######## 
def setParamMatrixDT(newVal) :
	global matrixDT

	if newVal == -1 :
		matrixDT = matrixDT
	elif newVal == -2 :
		matrixDT = np.float32
	else :
		matrixDT = newVal
	#end if

	return
#end def ######## ######## ######## 





######## ######## ######## ########
# ERROR CHECK: verify directory exists
# Input ----
#   path, str: path to save the file
#	create, boolean: whether to create missing dir
#	quiet, boolean: whether to quietly return T/F
# Returns ----
#	exists, boolean: indicates existence of directory
def verifyDirectory(path, create, quiet) :
	exists = True
	if not os.path.isdir(path) :
		if create :
			print("Creating path: {}".format(path))
			os.makedirs(path)
		elif quiet :
			exists = False
		else :
			print ( "ERROR: Specified path doesn't exist:" +
				" {}".format(path) )
			sys.exit()
	#end if
	return exists
#end def ######## ######## ######## 
######## ######## ######## ########
# ERROR CHECK: verify file exists
# Input ----
#   path, str: path to save the file
#	name, str: name of the file (w/ extension)
#	create, boolean: whether to create missing dir
#	quiet, boolean: whether to quietly return T/F
# Returns ----
#	exists, boolean: indicates existence of file
def verifyFile(path, name, quiet) :
	exists = True

	# First check the directory
	if not (path == '') :
		exists = verifyDirectory(path, False, quiet)

	# Then look for the file
	if not path.endswith('/') :
		path = path+'/'
	#end if
	if not os.path.isfile(path+name) :
		if quiet:
			exists = False
		else :
			print ( "ERROR: Specified file doesn't exist:" +
				" {}".format(path) )
			sys.exit()
	#end if

	return exists
#end def ######## ######## ######## 





######## ######## ######## ########
# Function: save a list to a text file
# Input:
#   path, str: path to save the file
#   name, str: name of file to save
#   theList, list of str - list of items to save
#       ASSUMPTION: list is already properly ordered
# Returns: nothing
# Creates: file containing ordered list of items
def saveListToText(path, name, theList) :

	# If folder doesn't exist, create it
	if not os.path.exists(path) :
		os.makedirs(path)
	#end if

	theFile = open(path + name, 'w')
	firstline = True
	for item in theList :
		if firstline :
			firstline = False
		else :
			theFile.write("\n")
		#end if
		theFile.write("{}".format(item))
	#end if
	theFile.close()

	return
#end def ######## ######## ######## 



######## ######## ######## ########
# Function: save the given matrix as a .txt file
# Input ----
#	matrix, (NxN) list: the values to save
#	mname, str: name of the file to save
#	mpath, str: path to the folder to save the file
#	integer, bool: True means save values as int()
# Returns ----
#	nothing
def saveMatrixText(matrix, mname, mpath, integer) :

	# If folder doesn't exist, create it
	if not os.path.exists(mpath) :
		os.makedirs(mpath)
	#end if

	# Open the file
	if mname.endswith('.txt') :
		fout = open(mpath + mname, "w")
	else :
		fout = open(mpath + mname + ".txt", "w")


	# Write to the file
	firstR = True
	for i in range(0, matrix.shape[0]) :

		# if not the first row, start with \n
		if firstR :
			firstR = False
		else :
			fout.write("\n")
		#end if

		firstC = True
		for j in range(0, matrix.shape[1]) :

			# if not the first col, start with \t
			if firstC :
				firstC = False
			else :
				fout.write("\t")
			#end if

			# Write the value to file
			#	If integer = True, write as an integer
			if integer :
				fout.write("{}".format( float(matrix[i,j]) ))
			else :
				fout.write("{}".format( matrix[i,j] ))
	#end loop
	fout.close()

	return
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Take the list of genes from the sample, return
#   a list of genes in the network, and left-out genes
# Input ----
#   path, str: path to the network files
#   name, str: name of the network to use
#   geneList, str list: genes in the sample
# Returns ----
#   inList & outList, str list: which genes from the
#       given list are found in network & which aren't
def checkGenesInNetwork(path, name, geneList) :

	fname = path + name + '_MetaPaths/genes.txt'

	# ERROR CHECK: verify file exists
	if not os.path.isfile(fname) :
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	# Read in the genes from the file
	geneSet = set()
	fg = open(fname, 'r')
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
# Function: Given a set of sample genes, 
# Pre-reqs: checkGenesInNetwork(path, name, geneList)
#           saveListToText(path, name, theList)
# Input ----
#   nPath, str: path to the network files
#   nName, str: name of the network to use
#   oPath, str: directory to create output
#   sample, str list: all genes in the sample
#   percent, int/float: percent of sample to conceal
# Returns ----
#   kGenes, str list: genes to be used as "known" test set
#   hGenes, str list: genes to be used as concealed set
# Creates ----
#   outPath/known.txt: list of known genes
#   outPath/concealed.txt: list of hidden genes
#   outPath/ignored.txt: list of genes not in network
def partitionSample(nPath, nName, oPath, sample, percent) :

	# Convert percent to fraction, if needed
	if (percent > 0) and (percent < 1) :
		# Assume percent is given as float (0, 1)
		pHide = percent
	elif (percent < 100) :
		# Assume percent is given as int [1, 100)
		pHide = percent / float(100)
	else :
		print ("ERROR: Invalid value given as the percent" +
			" of sample to hide.")
	#end if

	# Remove any genes not in the network
	inGenes, outGenes = checkGenesInNetwork(nPath, nName, sample)

	# Number of genes to keep vs hide
	numHide = int(len(inGenes) * pHide)
	numKeep = len(inGenes) - numHide

	# Separate the known and hidden genes
	random.shuffle(inGenes)
	kGenes = inGenes[0:numKeep]
	hGenes = inGenes[numKeep:len(inGenes)]

	# Save the gene lists to files
	kGenes.sort()
	saveListToText(oPath, fnKnownGenes, kGenes)
	hGenes.sort()
	saveListToText(oPath, fnConcealedGenes, hGenes)
	saveListToText(oPath, fnLeftOutGenes, outGenes)

	return kGenes, hGenes
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Read in the key.txt file regarding the 
#   metapath matrices
# Input ----
#   path, str: path to the network files
#   name, str: name of the network to use
# Returns ----
#   keyDict, dict
#       key, str: name of metapath
#       value, tuple: int is matrix/file ID number
#           bool where True means use matrix transpose
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
	fk = open(fname, "r")
	firstline = True
	for line in fk :

		# skip the first line
		if firstline :
			firstline = False
			continue
		#end if

		# separate the values
		line = line.rstrip()
#		print(line)
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
# Function: Read in the gene file. File is an ordered list
#   of genes, where the row number (starting at zero)
#   corresonds to the index in the matrix/list/etc where
#   the gene can be found.
# Input ----
#   fname, str: path & name to keep file
# Returns ----
#   iDict, dict: 
#       key, str: gene name as read from file
#       value, int: index to corresponding array
def readFileAsIndexDict(fname) :

	# ERROR CHECK: verify file exists
	if not os.path.isfile(fname) :
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	# Build the dictionary from the text file
	iDict = dict()
	gf = open(fname, "r")
	index = 0
	for line in gf :
		gene = line.rstrip()    # remove "\n"
		iDict[gene] = int(index)
		index += 1
	#end loop

	return iDict
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Read in the genes.txt file containing the 
#   gene-name headers to the metapath matrices
# Pre-reqs: readFileAsIndexDict(fname)
# Input ----
#   path, str: path to the network files
#   name, str: name of the network to use
# Returns ----
#   gDict, dict
#       key, str: name of gene
#       value, int: row/col index for that gene
def readGenesFile(path, name) :

	fname = path + name + "_MetaPaths/genes.txt"

	# The item to return
	gDict = readFileAsIndexDict(fname)

	return gDict
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Take the list of genes from the sample, return
#   a list of genes in the network, and left-out genes
# Input ----
#   numRows, int: number of random samples to create
#   numCols, int: number of items in each random sample
#   numItems, int: number of items in the source
#       (ie: the list from which the rand sample is pulled)
# Returns ----
#   rSamples, int array: 
def createRandomSamplesArray(numRows, numCols, numItems) :

	# The item to return
	rSamples = np.empty([numRows, numCols], dtype=np.int32)

	# Choose the random samples
	for r in range(0, numRows) :
		rSamples[r,:] = random.sample(range(numItems),
			numCols)
	#end loop

	return rSamples
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Read in a file as a line-by-line list of items
# Input ----
#   fname, str: path + name of the the sample files
# Returns ----
#   fItems, str list: ordered list of items from file
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
	fn = open(fname, "r")
	for line in fn :
		fItems.append( line.rstrip() )
	#end loop
	fn.close()

	fItems.sort()
	return fItems
#end def ######## ######## ######## 




######## ######## ######## ######## 
# Function: Read in the dataset from a samplename
#   Check for variants: ".txt", "_UP.txt", "_DN.txt"
# Pre-reqs: readFileAsList(fname)
# Input ----
#   fname, str: path + name of the the sample files
#   up, down: boolean -- only read the _UP or _DN file if true
# Returns ----
#   sNodes, str list: ordered list of names from file(s)
def readSampleFiles(sfile, up, down) :

	# The list of items to return
	sNodes = list()

	# Flag indicates a file existed and was read
	exists = False

	# First look for the file as named (no _UP or _DN)
	if os.path.isfile(sfile + ".txt") :
		temp = readFileAsList(sfile + ".txt")
		sNodes.extend(temp)
		exists = True
	#end if

	# Look for the _DN file
	if down and os.path.isfile(sfile + "_DN.txt") :
		temp = readFileAsList(sfile + "_DN.txt")
		sNodes.extend(temp)
		exists = True
	#end if

	if up and os.path.isfile(sfile + "_UP.txt") :
		temp = readFileAsList(sfile + "_UP.txt")
		sNodes.extend(temp)
		exists = True
	#end if

	# Alert user if nothing was read in
	if not exists :
		print("WARNING: no file found: {}".format(sfile))

	# Do NOT return duplicates
	uNodes = np.unique(sNodes) # sorted list of unique items
	return uNodes
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Given a list of names and appropriate dict,
#   convert to corresponding index values
# Input ----
#   names, str list: names of items to convert
#   iDict, dict
#       key, str: names of items
#       value, int: corresponding index values
# Returns ----
#   indices, int list: index values of the input items
def convertToIndices(names, iDict) :

	# The item to return
	indices = list()

	for name in names :
		indices.append(iDict[name])
	#end loop

	return indices
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Load the matrix containing the number of paths
#   of this type which join the nodes in the network
# Input ----
#   mpTuple [int, bool]: indicates which matrix file to use
#   path, str: path to the network files
#   name, str: name of the network to use
# Returns ----
#   matrix, int array: num paths between node pairs
def getPathMatrix(mpTuple, path, name, sizeOf) :

	zpad = fnMatrixZPad
#   fname = (path + name + "_MetaPaths/" +
#       "{}.npy".format(str(mpTuple[0]).zfill(zpad)) )

	prename = (path + name + "_MetaPaths/" +
		"{}".format(str(mpTuple[0]).zfill(zpad)) )
	if os.path.isfile(prename + '.gz') :
		fname = (path + name + "_MetaPaths/" +
		"{}.gz".format(str(mpTuple[0]).zfill(zpad)) )
	elif os.path.isfile(prename + '.txt') :
		fname = (path + name + "_MetaPaths/" +
		"{}.txt".format(str(mpTuple[0]).zfill(zpad)) )
	else :
		# ERROR CHECK: verify file exists
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	# Load the matrix
#   matrix = np.load(fname)
#	matrix = np.loadtxt(fname)

	# Declare the matrix
	if speedVsMemory :
		matrix = np.zeros([sizeOf, sizeOf])
	else :
		matrix = np.zeros([sizeOf, sizeOf], dtype=matrixDT)
	#end if

	# Read in the file, placing values into matrix
	row = 0
	with gzip.open(fname, 'rb') as fin :
		for line in fin :
			line = line.rstrip()
			ml = line.split()
			matrix[row,:] = ml[:]
			row += 1
	#end with

	# Convert to transpose if flag==True
	if mpTuple[1] :
		return np.transpose(matrix)
	else :
		return matrix
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Load the matrix containing the PathSim
#   metric for this meta-path
# Input ----
#   mpTuple [int, bool]: indicates which matrix file to use
#   path, str: path to the network files
#   name, str: name of the network to use
#	sizeOf, int: matrix dimensions (square)
#		NOTE: if sizeOf == -1, will find the dims
# Returns ----
#   matrix, int array: num paths between node pairs
def getSimMatrix(mpTuple, path, name, sizeOf) :

	zpad = fnMatrixZPad
#   fname = (path + name + "_MetaPaths/" +
#       "{}.npy".format(str(mpTuple[0]).zfill(zpad)) )

	prename = (path + name + "_MetaPaths/" +
		"Sxy-{}".format(str(mpTuple[0]).zfill(zpad)) )
	if os.path.isfile(prename + '.gz') :
		fname = prename + '.gz'
#		fname = (path + name + "_MetaPaths/" +
#		"Sxy-{}.gz".format(str(mpTuple[0]).zfill(zpad)) )
	elif os.path.isfile(prename + '.txt') :
		fname = prename + '.txt'
#		fname = (path + name + "_MetaPaths/" +
#		"Sxy-{}.txt".format(str(mpTuple[0]).zfill(zpad)) )
	else :
		# ERROR CHECK: verify file exists
		print ( "ERROR: Specified file doesn't exist:" +
			" {}.(gz/txt)".format(prename) )
		sys.exit()
	#end if

	# Load the matrix
#   matrix = np.load(fname)
#	matrix = np.loadtxt(fname)

	if sizeOf == -1 :
		# Get expected matrix size
	#TODO: pack this into a function
		fname = (path + name + "_MetaPaths/" +
			"{}.gz".format(str(0).zfill(fnMatrixZPad)) )
		sizeOf = 0
		with gzip.open(fname, 'rb') as fin :
			for line in fin :
				sizeOf += 1
		#end with
	#end if


	# Declare the matrix
	if speedVsMemory :
		matrix = np.zeros([sizeOf, sizeOf])
	else :
		matrix = np.zeros([sizeOf, sizeOf], dtype=matrixDT)
	#end if

	# Read in the file, placing values into matrix
	row = 0
	with gzip.open(fname, 'rb') as fin :
		for line in fin :
			line = line.rstrip()
			ml = line.split()
			matrix[row,:] = ml[:]
			row += 1
	#end with

	# Convert to transpose if flag==True
	if mpTuple[1] :
		return np.transpose(matrix)
	else :
		return matrix
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Find the number of paths of this type joining
#   the nodes in the sample
# Input ----
#   mpDict, {str: [int, bool]} dict:
#       key, str - name of the metapath
#       value, [int, bool] - which matrix file to use, and 
#           whether to use the transpose (inverse path)
# Returns ----
#   mpList, str list: ordered names of paths available,
#       less paths that are mirror-images of another
def removeInvertedPaths(mpDict) :

	# The item to return
	mpList = list()

	# Check the keys in the dict
	for key in mpDict.keys() :
		# If the boolean is True, then the path is an
		#   inverse of another; only append if false
		if mpDict[key][1]==False :
			mpList.append(key)
	#end loop

	mpList.sort()
	return mpList
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Find the number of paths of this type joining
#   the nodes in the sample
# Input ----
#   sample, int list: indices of the nodes in the sample
#   matrix, int array: num paths between node pairs
# Returns ----
#   count, int: num paths of this type in this sample
def getPathCountOne(sample, matrix) :

#TODO: If there is no difference in counts between a path
#   and it's inverse, then the inverse path need not be
#   calculated.
#NOTE: If a particular path is chosen for prediction, then
#   its inverse should be evaluated as well. (simulaneously)

	# Count the paths in the sample
	count = 0
#	for a in sample :
#		for b in sample :
#			if a != b :     # skip if a==b
#				count += matrix[a,b]
#	#end loop

	tempRows = matrix[sample,:]
	tempMatrix = tempRows[:,sample]
	count = tempMatrix.sum()
#	count = tempMatrix.sum() - np.sum(tempMatrix.diagonal())
#TODO: Should I remove links from node to self?

	return count
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Find the number of paths of this type joining
#   the nodes in the sample
# Input ----
#   samples, int 2D array: indices of the nodes in the sample
#       rows: individual samples
#   matrix, int array: num paths between node pairs
# Returns ----
#   counts, int list: num paths of this type in each sample
def getPathCountList(samples, matrix) :

#	counts = list()
	counts = np.zeros( (samples.shape[0]) )
	for i in range(samples.shape[0]) :
		# Count the paths in the sample
#		count = 0
#		for a in samples[i,:] :
#			for b in samples[i,:] :
#				if a != b :     # skip if a==b
#					count += matrix[a,b]
#		#end loop
#		counts.append(count)

		tempRows = matrix[samples[i,:],:]
		tempMatrix = tempRows[:,samples[i,:]]
		counts[i] = tempMatrix.sum()
#		counts[i] = tempMatrix.sum() - np.sum(tempMatrix.diagonal())
#TODO: Should I remove links from node to self?
	#end loop

	return counts
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: For a set of samples, find the mean and standard
#   devation of path counts for a specific path
# Input ----
#   rSamples, int 2D array: indices of the nodes in the sample
#       rows: individual samples
#   matrix, int array: num paths between node pairs
# Returns ----
#   mean, float list: mean count of paths for each sample
#   stdev, float list: standard deviation of counts per samp
def getPathMeans(rSamples, matrix) :

	# Get the path counts for each random sample
	counts = getPathCountList(rSamples, matrix)

	# Return the mean  and standard deviation
	mean = np.mean(counts)
	stdev = np.std(counts)
	return mean, stdev
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: For a set of samples and random counter-samples,
#	return the percent of random samples for which the
#	original had a higher count
#	(the percent of random samples beat by the original)
# Input ----
#   oCount, int: number of times this path exists in test sample
#   rSamples, int 2D array: indices of the nodes in the sample
#       rows: individual samples
#   matrix, int array: num paths between node pairs
# Returns ----
#   mean, float list: mean count of paths for each sample
#   stdev, float list: standard deviation of counts per samp
def getPercentile(oCount, rSamples, matrix) :

	# Get the path counts for each random sample
	rCounts = getPathCountList(rSamples, matrix)

	# Count how many times the original sample is more
	#   dense than a random sample
	numBeat = 0
	rankSum = 0
	for c in rCounts :
		rankSum += (oCount - c) / float(oCount + 0.001)
		if oCount > c :
			numBeat += 1
	#end loop
#TODO: test this:
#	numBeat = sum([1 for c in rCounts if (oCount > c)])
	percentile = numBeat / float(len(rSamples)) * 100
	rankAvg = rankSum / float(len(rCounts))

	return percentile, rankAvg
#end def ######## ######## ######## 



# ######## ######## ######## ######## 
# # Function: For a set of samples and random counter-samples,
# #	return the average percent difference between original
# #	and random. ie: if the randoms typically had 90% of the
# #	count as the original, return +0.1; if 110%, return -0.1
# # Input ----
# #   oCount, int: number of times this path exists in test sample
# #   rSamples, int 2D array: indices of the nodes in the sample
# #       rows: individual samples
# #   matrix, int array: num paths between node pairs
# # Returns ----
# #   mean, float list: mean count of paths for each sample
# #   stdev, float list: standard deviation of counts per samp
# def getPercentDifference(oCount, rSamples, matrix) :

# 	# Get the path counts for each random sample
# 	rCounts = getPathCountList(rSamples, matrix)

# 	# Count how many times the original sample is more
# 	#   dense than a random sample
# 	rankSum = 0
# 	for c in rCounts :
# 		rankSum += (oCount - c) / float(oCount + 0.001)
# 	#end loop
# 	rankAvg = rankSum / float(len(rCounts))

# 	return rankAvg
# #end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Calculate the statistics necessary to find
#   outstanding paths within the sample
# Pre-reqs: removeInvertedPaths(mpDict)
#           getPathMatrix(mpTuple, path, name)
#           getPathCountOne(sample, matrix)
#           getPercentile(oCount, rSamples, matrix)
#           getPathMeans(rSamples, matrix)
# Input ----
#   samples, int 2D array: indices of the nodes in the sample
#       rows: individual samples
#   mpDict, {str: [int, bool]} dict:
#       key, str - name of the metapath
#       value, [int, bool] - which matrix file to use, and 
#           whether to use the transpose (inverse path)
#   path, str: path to the network files
#   name, str: name of the network to use
# Returns ----
#   sCount, int list: num paths of this type in this sample
#   rMeans, float list: mean count of path in rand samples
#   rStDev, float list: standard deviation of rand samples
#   zScore, float list: z-Score for each path type
def calculateStatistics(sample, rSamples, mpDict,
	path, name) :

	# The items to return
	sCount = list()     # num of each path in given sample
	rMeans = list()     # mean count of e. p. in rand samples
	rStDev = list()     # stand dev of e. p. in rand samples
	zScore = list()     # z-Score of e. p. in rand samples
	percents = list()   # percentile rank of samp for e. p.

	# An iterable list of metapaths
	mpList = removeInvertedPaths(mpDict)

	# Get expected matrix size
#TODO: pack this into a function
	fname = (path + name + "_MetaPaths/" +
		"{}.gz".format(str(0).zfill(fnMatrixZPad)) )
	sizeOf = 0
	with gzip.open(fname, 'rb') as fin :
		for line in fin :
			sizeOf += 1
	#end with

	count = int(1)
	# Calculate the stats for each metapath
	for mp in mpList :

		matrix = getPathMatrix(mpDict[mp], path, name, sizeOf)

		tCount = getPathCountOne(sample, matrix)
		sCount.append(tCount)

		tPercent = getPercentile(tCount, rSamples, matrix)
		percents.append(tPercent)

		tMeans, tStDev = getPathMeans(rSamples, matrix)
		rMeans.append(tMeans)
		rStDev.append(tStDev)

		tScore = (tCount - tMeans) / (tStDev + 0.0001)
		zScore.append(tScore)

		if verbose :
			if not (count % 15) :
				print ("    tested {} paths".format(count))
		#end if
		count += 1
	#end loop

	return sCount, rMeans, rStDev, zScore, percents
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: choose an unused name for the output path
# Input ----
#   path, str: path to where output should be saved
#   dirPre, str: prefix of the folder name to return
# Returns ----
#   dirFull, str: name of output file (without path)
def nameOutputPath(path, dirPre) :

#TODO: replace with the verify function
	# ERROR CHECK: verify directory exists
	if not os.path.isdir(path) :
		print ( "ERROR: Specified path doesn't exist:" +
			" {}".format(path) )
		sys.exit()
	#end if

	zpad = fnOutputZPad

	# Set of all sub-folders in the path
	dirSet = [name for name in os.listdir(path)
		if os.path.isdir(path + name)]

	# increment folder name until an unused one is found
	num = int(0)
	dirFull = dirPre + "-{}".format(str(num).zfill(zpad))
	while dirFull in dirSet :
		num += 1
		dirFull = dirPre + "-{}".format(str(num).zfill(zpad))
	#end loop
	dirFull = dirFull + '/'

	return dirFull
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: choose an unused name for the output file
# Input ----
#   path, str: path to the network files
#   name, str: name of the network to use
# Returns ----
#   fname, str: name of output file (without path)
def nameOutputFile(path, name) :

#TODO: replace with the verify function
	# ERROR CHECK: verify directory exists
	if not os.path.isdir(path) :
		print ( "ERROR: Specified path doesn't exist:" +
			" {}".format(path) )
		sys.exit()
	#end if

	zpad = fnOutputZPad

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
# Pre-reqs: removeInvertedPaths(mpDict)
#           nameOutputFile(path, name)
# Input ----
#   path, str: path to the output files
#   name, str: name to use in creating output file
#   nName, str: name of the network used
#   sName, str: name of the sample used
#   mpDict, {str: [int, bool]} dict:
#       key, str - name of the metapath
#       value, [int, bool] - which matrix file to use, and 
#           whether to use the transpose (inverse path)
#   counts, int list: num paths of this type in this sample
#   means, float list: mean count of path in rand samples
#   stdevs, float list: standard deviation of rand samples
#   scores, float list: z-Score for each path type
#   leftout, str list: list of genes in original sample
#       which don't appear in the network
# Returns ----
#   fname, str: name of output file (without path)
def writeOutputOneSample(path, name, nName, sName, mpDict, 
	counts, means, stdevs, scores, percents, leftout) :

	delim = textDelim

	# Name the output file
	fname = nameOutputFile(path, name)

	# Get the ordered list of metapaths
	mpList = removeInvertedPaths(mpDict)

	# Write the file header
	fn = open(path+fname, "wb")
	fn.write("The occurrence of metapaths with in a sample ...\n")
	fn.write("network:{}{}\n".format(delim, nName))
	fn.write("sample:{}{}\n".format(delim, sName))
	fn.write("\n")

	# Write the body of the file
	fn.write("sample{}random samples\n".format(delim))
	fn.write("count{0}mean{0}std dev{0}".format(delim) +
		"z-score{0}percentile{0}path name{0}length\n".format(delim))
	i = 0
	for mp in mpList :
		length = mp.count('-') + 1
		fn.write("{1:d}{0}{2}{0}{3}{0}{4}{0}{5}{0}{6}{0}{7}\n".format(
			delim, int(counts[i]), means[i], stdevs[i],
			scores[i], percents[i], mp, length) )
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



######## ######## ######## ######## 
# Function: Use a chosen statistic to select the top K
#   metapaths for use in making predictions
# Input ----
#   k, int: number of paths to return
#   ranker, float list: stat ordered by mp name
#   mpDict, {str: [int, bool]} dict:
#       key, str - name of the metapath
#       value, [int, bool] - which matrix file to use, and 
#           whether to use the transpose (inverse path)
# Returns ----
#   kPaths, str list: the recommended top K paths
def chooseTopKPathsSimple(k, ranker, mpDict) :

	# The ordered list of metapaths
	mpList = removeInvertedPaths(mpDict)

	rankMax = np.amax(ranker)

	# Create an array of the paths & stats to rank
	pathArray = np.zeros( len(ranker), dtype=[('pIndex', 'i4'),
		('stat', 'f2'), ('pLength', 'u1'), ('stat_inverse', 'f2')] )

	for i in range(len(ranker)) :
		length = mpList[i].count('-') + 1
		invertRank = rankMax - ranker[i]
		pathArray[i] = (i, ranker[i], length, invertRank)
	#end loop

#TODO: sort by ranker DESCENDING and by length ASCENDING
# Can I do this, or do I need to invert ranker?

	# Sort array DESCENDING by rank and ASCENDING by length
	pathArray = np.sort(pathArray, order=['stat_inverse', 'pLength'])

	# Convert to lists to return
	topPaths = [mpList[idx] for idx in pathArray['pIndex'][0:k]]
	topRanks = [value for value in pathArray['stat'][0:k]]

	return topPaths, topRanks
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Calculate the Group PathSim for all genes not
#	in the test sample
# Input ----
#	path, str: path to the network files
#	name, str: name of the network to use
#	mpTuple [int, bool]: indicates which matrix file to use
#	sample, int list: indices of the genes in the sample
# Returns ----
#	scores, float list: G P S for each gene not in sample
def applyGroupPathSim(path, name, mpTuple, sample) :

	# Get expected matrix size
#TODO: pack this into a function
	fname = (path + name + "_MetaPaths/" +
		"{}.gz".format(str(0).zfill(fnMatrixZPad)) )
	sizeOf = 0
	with gzip.open(fname, 'rb') as fin :
		for line in fin :
			sizeOf += 1
	#end with

	# Get the corresponding path matrix
	matrix = getPathMatrix(mpTuple, path, name, sizeOf)

	# Apply the modified Group PathSim

	# Number of paths connecting sample to self:
	temp1 = matrix[sample,:]
	temp2 = temp1[:,sample]
	pXX = np.sum(temp2)
	del temp1, temp2

	# The item to return
	scores = list()

	# Get the indices for all genes not in the test sample
	notSample = ([n for n in range(matrix.shape[0])
		if n not in sample])
	notSample.sort()	# maybe redundant, but want to be sure

	for i in notSample :
		# Number of paths connecting gene to self:
		pYY = matrix[i,i]
		# Number of paths from sample to gene:
		pXY = np.sum(matrix[sample,i])
		# Number of paths from gene to sample:
		pYX = np.sum(matrix[i,sample])

# TODO: If (pXX + pYY) == 0 ... set = 1 ?
		scores.append( (pXY + pYX) / (pXX + pYY + 0.0001) )
	#end if

#	# For intermediate verification:
#	geneList = geneIndex.keys()
#	geneList.sort()
#	print geneList
#	print "  name  index   score"
#	for i in range(len(sampOutdex)) :
#		print "  {}, {}  :  {}".format(geneList[sampOutdex[i]],
#			sampOutdex[i], geneRanks[i])
#	#end if

	return scores
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Write the file listing the items & rankings
# Input ----
#	path, str: directory to write output file
#	statArray, float array: the stat used to rank the items
#	itemDict, {str: int} dict:
#		key, str - name of the item
#		value, int - index of the item within a sorted list
#	itemIndex, int list: indices of items (known/test) against
#		which the rest of the items (hidden/concealed) are ranked
# Returns ----
#	filename, str: path & filename to output file
# Creates ----
#	pathlist-###.txt: original version of the output file
def writeItemRanks(path, statArray, itemDict, itemIndex, colHeader) :


	statSum = np.sum(statArray, axis=1)
	statAvg = np.mean(statArray, axis=1)

	# Write the file header
	delim = '\t'
	filename = path + 'scores.txt'
	fs = open(filename, 'wb')
	fs.write('')
	fs.write('{}Genes scored by similarty to sample ...\n'.format(textDelim))
#	fs.write('network:{}{}\n'.format(delim, eName))
#	fs.write('sample:{}{}\n'.format(delim, sName))
#	fs.write('known:{}{}\n'.format(delim, len(gKnown)))
#	fs.write('concealed:{}{}\n'.format(delim, len(gHidden)))
	fs.write('\n')


	# Get the indices for all genes not in the test sample
	itemOutdex = [n for n in range(len(itemDict)) if n not in itemIndex]
	#print itemOutdex
	geneList = itemDict.keys()
	geneList.sort()
	geneList = [geneList[g] for g in itemOutdex]
	del itemOutdex
	#print geneList


	#fs.write('{}'.format(delim))
	for bp in colHeader :
		fs.write('{}{}'.format(delim, bp))
	#end loop
	fs.write('{}sum{}avg\n'.format(delim, delim))

	#print statArray.shape, len(statSum), len(statAvg)

	for i in range( statArray.shape[0] ) :
		fs.write('{}'.format(delim))
		for j in range( statArray.shape[1] ) :
			fs.write( '{}{}'.format(statArray[i,j], delim) )
		#end loop
	#	print i
		fs.write( '{}{}{}{}{}\n'.format(statSum[i], delim,
			statAvg[i], delim, geneList[i]) )
	#end loop
	fs.write('\n')

	return filename
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Choose N random samples, return as an array of
#	gene indices, ! Use node-binning !
# Input ----
#	path, str: path to the network files
#	name, str: name of the network to use
#	sample, str list: list of nodes by name
#		ASSUME: only contains nodes in the network
#	indexDict, dict:
#       key, str: gene name as read from file
#       value, int: index to corresponding array
#	cutoffPercent, float list: values to use as bin cutoffs
#		ex: [0.4] means lower bin contains 40% of lowest-degree nodes
#	eType, str: the desired edge type along which to bin,
#		'all' means total degree count
#	numRandSamples, int: desired #, indicates # rows
# Returns ----
#   rSamples, int array: 
def createRandomSamplesBinned(path, name, sample,
	indexDict, cutoffPercent, eType, numRandSamples) :

#	print "verbosity = {}".format(verbose)

	fname = path + name + '/node-degree.txt'

	# ERROR CHECK: verify file exists
	if not os.path.isfile(fname) :
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	if verbose :
		print ("Opening {}".format(fname))


	# Read in the file:
	dfile = open(fname, 'r')
	col = -1
	degreeDict = dict()
	for line in dfile :
		line = line.rstrip()
		lv = line.split(textDelim)

		# Locate the column corresponding to eType
		if lv[0] == 'HEADER' :
			for i in range(len(lv)) :
				if lv[i] == eType :
					col = i
			# end loop
			if verbose :
				print ("    binning along edge type: {} (column {})".format(eType, col))
			if col == -1 :
				print ( "ERROR: Specified edge type doesn't appear" +
					" in node-degree.txt" )
				sys.exit()
			#end if

		# Once column is found, create the node-degree dict
		else :
			degreeDict[lv[0]] = int(lv[col])
	#end loop

#	print(degreeDict)

	# Define the bin cutoffs
	values = list(degreeDict.values())
	values.sort()
	cutoffs = list()
	for cp in cutoffPercent :

		if (cp <= 1) and (cp >= 0) :
			decimal = cp
		elif cp < 100 :
			decimal = cp / 100
		else :
			print ( "ERROR: Specified node-binning cutoff percent" +
				" is outside range: {}".format(decimal) )
		#end if

		# Find the value of the gene at cp %
		# Goal: capture X % of nodes, sorted by degree
		location = int(round(cp * len(values)))
#		print("{} = {}*{}".format(location, cp, len(values)))
		cutoffs.append(values[location])
	#end loop
	cutoffs.append(np.amax(values) + 1)	# The 100% ceiling


	# Bin the sample
	distribution = [0] * (len(cutoffs) + 1)
	for node in sample :
		# Get the sample node's degree
		this = degreeDict[node]

		# Increment the appropriate bin
		for i in range(len(cutoffs)) :
			if this < cutoffs[i] :
				distribution[i] += 1
				break
		# end loop
	#end loop

	if verbose :
		print ("    degree distribution of sample: {}".format(distribution))
		print ("    with cutoffs: {}".format(cutoffs))
	#end if



	# Build inverse array as {bin: [genes]}
	binDict = dict()
	for i in range(len(distribution)) :
		binDict[i] = list()
	#end loop

	# Place the node names into the bins
	for key in degreeDict :
		# Add key to lowest matching bin number
		for i in range(len(cutoffs)) :
			if degreeDict[key] < cutoffs[i] :
				binDict[i].append(key)
				break
		#end loop
	#end loop
	del degreeDict	# Your services are no longer needed.

	if verbose :
		distBins = list()
		for i in range(len(binDict)) :
			distBins.append(len(binDict[i]))
		print ("    degree distribution of network: {}".format(distBins))
	#end if


	# Populate the random samples matrix (rows = samples, cols = indices to nodes)
	rSamples = np.empty([numRandSamples, len(sample)], dtype=np.int32)
	for row in range(0, numRandSamples) :
		oneRSample = list()

		# Select genes from each bin, add to random sample
		for i in range(len(binDict)) :

			if distribution[i] > len(binDict[i]) :
				print ( "ERROR: For bin {} the number of nodes in sample ({})" +
					" exceeds number for the network {}".format(i, 
						distribution[i], len(binDict[i])) )
			#end if

			# get the indices from the list of nodes in bin
			indices = random.sample(range(len(binDict[i])), distribution[i])
			indices.sort()
			# convert to list of indices from the passed dict
			tempList = [indexDict[binDict[i][idx]] for idx in indices]

			oneRSample.extend(tempList)
		#end loop

		rSamples[row,:] = oneRSample
	#end loop


	return rSamples
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Write a simple ranked_genes.txt file
# Input ----
#	path, str: directory to write output file
#	statArray, float array: the stat used to rank the items
#	itemDict, {str: int} dict:
#		key, str - name of the item
#		value, int - index of the item within a sorted list
#	itemIndex, int list: indices of items (known/test) against
#		which the rest of the items (hidden/concealed) are ranked
#	hiddenSet, str list: names of the genes which were concealed
#	cutoffs, int list: number of genes which might be returned at once
# Returns ----
#	nothing
# Creates ----
#	ranked_genes.txt: original version of the output file
#	returned_cutoffs.txt: shows how many desired genes exist within N predicted
def writeRankedGenes(path, statArray, itemDict, itemIndex, hiddenSet, cutoffs) :

	rankedFile = 'ranked_genes.txt'
	cutoffFile = 'returned_cutoffs.txt'

	# Create the numpy record/structured array
	rankList = np.recarray(statArray.shape[0],
		dtype=[('inverse', 'f4'), ('score', 'f4'), ('names', nodeDT)])
	# The scores come from the average of the score array
	statAvg = np.mean(statArray, axis=1)
	rankList['score'] = statAvg
	rankList['inverse'] = 1 - statAvg
	# The gene names exclude the known/test sample
#	itemOutdex = [n for n in range(len(itemDict)) if not in itemIndex]
#	itemOutdex.sort()
	itemNames = list(itemDict.keys())
	itemNames.sort()
#	rankList['index'] = itemNames[itemOutdex]
#	rankList['names'] = itemNames[n for n in range(len(itemDict)) if n not in itemIndex]
	rankList['names'] = [itemNames[n] for n in range(len(itemDict)) if n not in itemIndex]

	# Sort by the score
	rankList.sort(order=['inverse', 'names'])


	# Open the files to write
	fouta = open(path+rankedFile, 'w')
	foutb = open(path+cutoffFile, 'w')

	# Write the header for the cutoffs file
	foutb.write("Number of True Positives returned in top N predicted...")
	foutb.write("\nReturned{}TruePos".format(textDelim))

	foundCount = 0
	firstLine = True
	for i in range(len(rankList)) :

		# Write the body of the ranked file
		if not firstLine :
			fouta.write("\n")
		firstLine = False
		fouta.write("{}{}{}".format(rankList['score'][i], textDelim, rankList['names'][i]))

		# Write the body of the cutoffs file
		if rankList['names'][i] in hiddenSet :
			foundCount += 1
		if i in cutoffs :
			foutb.write("\n{}{}{}".format(i, textDelim, foundCount))

	#end loop
	fouta.close()
	foutb.close()

	return
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Write a simple ranked_genes.txt file
# Input ----
#	path, str: directory to write output file
#	suffix, str: suffix to the filename
#	statArray, float array: the stat used to rank the items
#	itemDict, {str: int} dict:
#		key, str - name of the item
#		value, int - index of the item within a sorted list
#	itemIndex, int list: indices of items (known/test) against
#		which the rest of the items (hidden/concealed) are ranked
#	hiddenSet, str list: names of the genes which were concealed
#	cutoffs, int list: number of genes which might be returned at once
# Returns ----
#	nothing
# Creates ----
#	ranked_genes.txt: original version of the output file
#	returned_cutoffs.txt: shows how many desired genes exist within N predicted
def writeRankedGenes02(path, suffix, statArray, itemDict, itemIndex, hiddenSet, cutoffs) :

	rankedFile = 'ranked_genes-' + suffix + '.txt'
	cutoffFile = 'returned_cutoffs-' + suffix + '.txt'

	# Create the numpy record/structured array
	rankList = np.recarray(statArray.shape[0],
		dtype=[('inverse', 'f4'), ('score', 'f4'), ('names', nodeDT)])
	# The scores come from the average of the score array

#	print(statArray)

#TODO: Mean vs Sum ... any difference ??
	statAvg = np.mean(statArray, axis=1)
	rankList['score'] = statAvg
	rankList['inverse'] = 0 - statAvg

#	statSum = np.sum(statArray, axis=1)
#	rankList['score'] = statSum
#	rankList['inverse'] = 0 - statSum

	# The gene names exclude the known/test sample
#	itemOutdex = [n for n in range(len(itemDict)) if not in itemIndex]
#	itemOutdex.sort()
	itemNames = list(itemDict.keys())
	itemNames.sort()

#	print(itemNames[0])

#	rankList['index'] = itemNames[itemOutdex]
#	rankList['names'] = itemNames[n for n in range(len(itemDict)) if n not in itemIndex]
	rankList['names'] = [itemNames[n] for n in range(len(itemDict)) if n not in itemIndex]
#TODO: Why is this writing as a byte stream, not str ??

	# Sort by the score
	rankList.sort(order=['inverse', 'names'])


	# Open the files to write
	fouta = open(path+rankedFile, 'w')
	foutb = open(path+cutoffFile, 'w')

	# Write the header for the cutoffs file
	foutb.write("Number of True Positives returned in top N predicted...")
	foutb.write("out of {} concealed genes".format(len(hiddenSet)))
	foutb.write("\nReturned{}TruePos".format(textDelim))

#	print("hidden {}".format(hiddenSet))
#	print(rankList['names'])

	foundCount = 0
	firstLine = True
	for i in range(len(rankList)) :

		# Write the body of the ranked file
		if not firstLine :
			fouta.write("\n")
		firstLine = False
#		print(rankList['score'][i])
#		fouta.write("{}{}{}".format(rankList['score'][i], textDelim, rankList['names'][i]))
		fouta.write("{}{}{}".format(rankList['score'][i],
			textDelim, rankList['names'][i]))

		# Write the body of the cutoffs file
		if rankList['names'][i] in hiddenSet :
			foundCount += 1
		if i in cutoffs :
			foutb.write("\n{}{}{}".format(i, textDelim, foundCount))

	#end loop
	fouta.close()
	foutb.close()

	return rankList['names']
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Save the topK chosen paths w/ scores
# Input ----
#	path, str: directory to write output file
#	suffix, str: suffix to the filename
#	pathList, str list: list of the paths to save
#	scoreList, float list: list of the respective scores
# Returns ----
#	nothing
# Creates ----
#	ranked_paths-<suffix>.txt: list of the paths used
def writeChosenPaths(path, suffix, pathList, scoreList) :

	if len(pathList) != len(scoreList) :
		print("ERROR: pathList ({}) and scoreList ".format(len(pathList)) +
			"({}) are different lengths.".format(len(scoreList)))
	#end if

	# Open the output file
	outName = 'ranked_paths-' + suffix + '.txt'
	fout = open(path + outName, 'w')

	# Write the header
	fout.write("Top paths & scores as chosen by the" +
		" {} method".format(suffix))

	# Write the body
	for i in range(len(scoreList)) :
		fout.write("\n{}{}{}".format(scoreList[i],
			textDelim, pathList[i]))
	#end loop

	fout.close()

	return
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Write a simple ranked_paths.txt file
# Input ----
#	path, str: directory to write output file
#   ranker, float list: stat ordered by mp name
#   mpDict, {str: [int, bool]} dict:
#       key, str - name of the metapath
#       value, [int, bool] - which matrix file to use, and 
#           whether to use the transpose (inverse path)
# Returns ----
#	nothing
# Creates ----
#	ranked_paths.txt: original version of the output file
def writeRankedPaths(path, suffix, ranker, mpDict) :

#	rankedFile = 'ranked_paths.txt'
	if not suffix :
		rankedFile = 'ranked_paths.txt'
	else :
		rankedFile = 'ranked_paths-' + suffix + '.txt'

	# The ordered list of metapaths
	mpList = removeInvertedPaths(mpDict)

	if len(mpList) != len(ranker) :
		print ( "WARNING: The list of paths ({}) ".format(len(mpList)) +
			"and scores ({}) are unequal length.".format(len(mpList)) )
	#end if


	# Create an array of the paths & stats to rank
	# 	dtype=object allows me to store variable-sized str
	pathArray = np.recarray( len(ranker), dtype=[('name', object),
		('stat', 'f4'), ('length', 'u1'), ('stat_inverse', 'f4')] )

	rankMax = np.amax(ranker)
	for i in range(len(ranker)) :
		length = mpList[i].count('-') + 1
		invertRank = rankMax - ranker[i]
		pathArray[i] = (mpList[i], ranker[i], length, invertRank)
	#end loop

	# Sort array DESCENDING by rank and ASCENDING by length
	pathArray = np.sort(pathArray, order=['stat_inverse', 'length', 'name'])


	# Write to file
	fout = open(path+rankedFile, 'w')
	firstLine = True
	for i in range(len(pathArray)) :
		if not firstLine :
			fout.write("\n")
		firstLine = False
		fout.write("{1}{0}{2}{0}{3}".format(textDelim, pathArray['stat'][i],
			pathArray['name'][i], pathArray['length'][i]))
	#end loop
	fout.close()


	return
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Read in the ranked_paths file
# Input ----
#   path, str: path to file
#	name, str: name of the file
# Returns ----
#	np structured array :
#		(percent, path name, length)
def readRankedPaths(path, name) :

	if not path.endswith('/') :
		path = path + '/'

	# ERROR CHECK: verify file exists
#	verifyFile(path, 'ranked_paths.txt', False)
	verifyFile(path, name, False)


	# Get the number of lines in the file
#	fname = path + 'ranked_paths.txt'
	fname = path + name
	nRows = 0
	with open(fname, 'r') as fin :
		nRows = sum(1 for line in fin)
#		for line in fin :
#			nRows += 1
	#end with


	# Create an array of the paths & stats to rank
	# 	dtype=object allows me to store variable-sized str
	pathArray = np.recarray( nRows, dtype=[('name', object),
		('stat', 'f4'), ('stat_inverse', 'f4'), ('length', 'u1')] )


	# Read in the file
	with open(fname, 'r') as fin :
		i = 0
		for line in fin :
			line = line.rstrip()    # remove "\n"
			lv = line.split(textDelim)

			pathArray[i] = (lv[1],
				float(lv[0]), (0.0 - float(lv[0])), int(lv[2]))
			i += 1
	#end with

	return pathArray
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Create a list of samples contained in folder
# Input ----
#	sPath, str: path where samples stored
# Returns ----
#	sNames, str list: sorted list of sample names
def getSampleNamesFromFolder(path) :

	verifyDirectory(path, False, False)

	# Get list of all text files in folder
	fNames = [f for f in listdir(path) if f.endswith('.txt')]
	if verbose :
		print (fNames)

	# Identify & create list of sample names in folder
	sNames = list()
	for item in fNames :
		# Strip the extension and any "_UP" or "_DN"
		newItem = item[:-4]
		if newItem.endswith('_UP') or newItem.endswith('_DN') :
			newItem = newItem[:-3]
		#end if
		sNames.append(newItem)
	#end loop

	sNames = np.unique(sNames)	# also sorts list
	return sNames
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Write a text file where the columns are given as lists
# Input ----
#	path, str: directory to write output file
#	fname, str: name of the file to write
#	columnList, list of str lists:
#		each entry in columnList represents a column
#		where each entry is a string to write to the file
# Returns ----
#	nothing
# Creates ----
#	ranked_paths.txt: original version of the output file
def writeGenericLists(path, fname, columnList) :

	verifyDirectory(path, True, False)

# ASSUME: the contained lists are of equal length

	fout = open(path+fname, 'w')

	for i in range(len(columnList[0])) :
		fout.write("{}".format(columnList[0][i]))

		for j in range(1, len(columnList)) :
			fout.write("{}{}".format(textDelim, columnList[j][i]))
		#end loop

		fout.write("\n")
	#end if

	fout.close()

	return
#end def ######## ######## ######## 




######## ######## ######## ######## 
# Function: Return list of paths to folders in the 
#	given root directory
# Input ----
#	root, str: path where the folders reside
# Returns ----
#	subDirs, str list: sorted list of subdirectories
#		contains full path: root+subdir
def getSubDirectoryList(root) :

	verifyDirectory(root, False, False)

	if not root.endswith('/') :
		root = root+'/'
	subDirs = [(root+d+'/') for d in os.listdir(root) if os.path.isdir(root+d)]
	subDirs.sort()

	return subDirs
#end def ######## ######## ######## 




######## ######## ######## ######## 
# Function: Return expected height/width of matrix
#	Used for pre-allocating memory
# Input ----
#	path, str: path to the matrix
#	name, str: matrix file name (w/ extension)
# Returns ----
#	nRows, nCols, int: number of rows & columns
# ASSUMES: no header/footer -- matrix is the entire file
def getMatrixDimensions(path, name) :

#	verifyDirectory(path, False, False)
	verifyFile(path, name, False)

	nRows = 0
	nCols = -1

	if name.endswith('.gz') :
		with gzip.open(path+name, 'rb') as fin :
			for line in fin :
				nRows += 1
				temp = len( line.split(textDelim) )
				if temp == 0 :
					break
				#end if
				if nCols == -1 :
					nCols = temp
				else :
					nCols = min(temp, nCols)
		#end with
	elif name.endswith('.txt') :
		with open(path+name, 'rb') as fin :
			for line in fin :
				nRows += 1
				temp = len( line.split(textDelim) )
				if temp == 0 :
					break
				#end if
				if nCols == -1 :
					nCols = temp
				else :
					nCols = min(temp, nCols)
		#end with
	else :
		print("  ERROR: Expected either '.gz' or '.txt' as extension: {}".format(name))
		sys.exit()
	#end if

	return nRows, nCols
#end def ######## ######## ######## 




######## ######## ######## ######## 
# Function: Read contents of file into a matrix
# Input ----
#	path, str: path where the file resides
#	name, str: name of the file
# Returns ----
#	matrix, numpy matrixDT: the 2D matrix
def readFileAsMatrix(path, name) :

	# Declare the matrix (get size)
	nRows, nCols = getMatrixDimensions(path, name)
	matrix = np.zeros([nRows, nCols], dtype=matrixDT)

	# Read in the file
	if name.endswith('.gz') :
		with gzip.open(path+name, 'rb') as fin :
			count = 0
			for line in fin :
				line = line.rstrip('\n')
				lv = line.split(textDelim)
				matrix[count,:] = lv[:]
				count += 1
		#end with
	elif name.endswith('.txt') :
		with open(path+name, 'rb') as fin :
			for line in fin :
				line = line.rstrip('\n')
				lv = line.split(textDelim)
				matrix[count,:] = lv[:]
				count += 1
		#end with
#	else :
#		print("  ERROR: Expected either '.gz' or '.txt' as extension: {}".format(name))
#		sys.exit()
	#end if

	return matrix
#end def ######## ######## ######## 