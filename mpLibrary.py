# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# functions used by the Metapath Analysis scripts
#
# ---------------------------------------------------------

import os.path
from os import listdir
import sys
import numpy as np
import random



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

# Length to pad the matrix file names:
keyZPad = 5
# Length to pad the output file names:
fnZPad = 3
# Data delimiter to use in the output file:
textDelim = '\t'
# Whether to print non-error messages within these funcs
verbose = True

######## ######## ######## ######## 






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

	theFile = open(path + name, 'wb')
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
	fg = open(fname, 'rb')
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
	gf = open(fname, "rb")
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
		rSamples[r,:] = random.sample(xrange(numItems),
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
	fn = open(fname, "rb")
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
		print "WARNING: no file found: {}".format(sfile)

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
def getPathMatrix(mpTuple, path, name) :

	zpad = keyZPad
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
	matrix = np.loadtxt(fname)

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
	for a in sample :
		for b in sample :
			if a != b :     # skip if a==b
				count += matrix[a,b]
	#end loop

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

	counts = list()
	for i in range(0, samples.shape[0]) :
		# Count the paths in the sample
		count = 0
		for a in samples[i,:] :
			for b in samples[i,:] :
				if a != b :     # skip if a==b
					count += matrix[a,b]
		#end loop
		counts.append(count)
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
# Function: For a set of samples, find the mean and standard
#   devation of path counts for a specific path
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
	for c in rCounts :
		if oCount > c :
			numBeat += 1
	#end loop
	percentile = numBeat / float(len(rSamples)) * 100

	return percentile
#end def ######## ######## ######## 



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

	# Calculate the stats for each metapath
	for mp in mpList :

		matrix = getPathMatrix(mpDict[mp], path, name)

		tCount = getPathCountOne(sample, matrix)
		sCount.append(tCount)

		tPercent = getPercentile(tCount, rSamples, matrix)
		percents.append(tPercent)

		tMeans, tStDev = getPathMeans(rSamples, matrix)
		rMeans.append(tMeans)
		rStDev.append(tStDev)

		tScore = (tCount - tMeans) / (tStDev + 0.0001)
		zScore.append(tScore)
	#end loop

	return sCount, rMeans, rStDev, zScore, percents
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: choose an unused name for the output file
# Input ----
#   path, str: path to where output should be saved
#   dirPre, str: prefix of the folder name to return
# Returns ----
#   dirFull, str: name of output file (without path)
def nameOutputPath(path, dirPre) :

	# ERROR CHECK: verify directory exists
	if not os.path.isdir(path) :
		print ( "ERROR: Specified path doesn't exist:" +
			" {}".format(path) )
		sys.exit()
	#end if

	zpad = fnZPad

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
# Function: Load the matrix containing the number of paths
#	of this type which join the nodes in the network
# Input ----
#	mpTuple [int, bool]: indicates which matrix file to use
#	path, str: path to the network files
#	name, str: name of the network to use
# Returns ----
#	matrix, int array: num paths between node pairs
def getPathMatrix(mpTuple, path, name) :

	zpad = keyZPad
#	fname = (path + name + "_MetaPaths/" +
#		"{}.npy".format(str(mpTuple[0]).zfill(zpad)) )

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
#	matrix = np.load(fname)
	matrix = np.loadtxt(fname)

	# Convert to transpose if flag==True
	if mpTuple[1] :
		return np.transpose(matrix)
	else :
		return matrix
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

	# Get the corresponding path matrix
	matrix = getPathMatrix(mpTuple, path, name)

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
# Function: Calculate the Group PathSim for all genes not
#	in the test sample
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
def writeItemRanks(path, statArray, itemDict, itemIndex, colHeader) :


	statSum = np.sum(statArray, axis=1)
	statAvg = np.mean(statArray, axis=1)

	delim = '\t'
	filename = path + 'scores.txt'
	fs = open(filename, 'wb')
	fs.write('')
	fs.write('Genes scored by similarty to sample ...\n')
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

	fname = path + name + '/node-degree.txt'

	# ERROR CHECK: verify file exists
	if not os.path.isfile(fname) :
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	if verbose :
		print "Opening {}".format(fname)


	# Read in the file:
	dfile = open(fname, 'rb')
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
				print "    {} at column {}".format(eType, col)
			if col == -1 :
				print ( "ERROR: Specified edge type doesn't appear" +
					" in node-degree.txt" )
				sys.exit()
			#end if

		# Once column is found, create the node-degree dict
		else :
			degreeDict[lv[0]] = int(lv[col])
	#end loop


	# Define the bin cutoffs
	values = degreeDict.values()
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
		print "    degree distribution of sample: {}".format(distribution)
		print "    with cutoffs: {}".format(cutoffs)
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
		print "    degree distribution of network: {}".format(distBins)
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
			indices = random.sample(xrange(len(binDict[i])), distribution[i])
			indices.sort()
			# convert to list of indices from the passed dict
			tempList = [indexDict[binDict[i][idx]] for idx in indices]

			oneRSample.extend(tempList)
		#end loop

		rSamples[row,:] = oneRSample
	#end loop


	return rSamples
#end def ######## ######## ######## 