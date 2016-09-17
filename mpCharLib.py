# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# functions used by the Metapath Characterization scripts
#
#
# Functions provided:
#	concatenatePaths(root, subDir)
#	
# ---------------------------------------------------------

import os.path
from os import listdir
import sys
import numpy as np
import gzip
#import random
#from sklearn import cluster as skc
#import warnings




######## ######## ######## ######## 
# PARAMETERS


######## ######## ######## ######## 




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
			os.makedirs(path)
			if not quiet :
				print("Creating path: {}".format(path))
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
	if path != '' :
		exists = verifyDirectory(path, False, quiet)
		# Then concatenate path & name
		if not path.endswith('/') :
			path = path + '/'
		fname = path + name
	else :
		fname = name
	#end if

	# Then look for the file
	if not os.path.isfile(fname) :
		if quiet:
			exists = False
		else :
			print ( "ERROR: Specified file doesn't exist:" +
				" {}".format(fname) )
			sys.exit()
	#end if

	return exists
#end def ######## ######## ######## 




######## ######## ######## ######## 
# Combine a root directory with a sub-directory
# Input ----
#	root, str: the root directory
#	subDir, str: the sub-directory
# Returns ----
#	path, str: the full combined path
def concatenatePaths(root, subDir) :
	if not root.endswith('/') :
		root = root + '/'
	if not subDir.endswith('/') :
		subDir = subDir + '/'
	#end if
	path = root + subDir

	return path
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Choose an unused name for the output path
# Input ----
#   path, str: path to where output should be saved
#   dirPre, str: prefix of the folder name to return
# Returns ----
#   dirFull, str: name of output file (without path)
def nameOutputPath(path, dirPre) :

	# ERROR CHECK: verify directory exists
	verifyDirectory(path, False, False)

	zpad = 3

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
			fout.write("{}{}".format('\t', columnList[j][i]))
		#end loop

		fout.write("\n")
	#end if

	fout.close()

	return
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

	verifyFile('', fname, False)

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
def getGeneDictionary(path, name) :

	fname = concatenatePaths(path, name)
	fname = fname + "genes.txt"

	# The item to return
	gDict = readFileAsIndexDict(fname)

	return gDict
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
def getPathDictionary(path, name) :

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
# Function: Load the matrix containing the number of paths
#   of this type which join the nodes in the network
# Input ----
#   mpTuple [int, bool]: indicates which matrix file to use
#   path, str: path to the network files
#   name, str: name of the network to use
# Returns ----
#   matrix, int array: num paths between node pairs
def getPathMatrix(mpTuple, path, name, sizeOf) :

#TODO: change name(s) from "get..." to "load..."

	speedVsMemory = False	# True favors speed, disables dtype
	matrixDT = np.float32	#TODO: any considerations here?
	fnMatrixZPad = 5


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
# Count number of rows in a metapath matrix
#	(num rows = num cols)
# Input ----
#   ePath, str: path to network
#	eName, str: folder containing processed network files
# Returns ----
#	mxSize, int: number of rows/columns in path matrix
def getPathMatrixSize(ePath, eName) :

	# number of figures in the zero-padded matrix file name
	fnMatrixZPad = 5

	# the item to return
	mxSize = 0

	# open and read through the file
	fname = (ePath+eName+"_MetaPaths/" +
		"{}.gz".format(str(0).zfill(fnMatrixZPad)) )
	with gzip.open(fname, 'rb') as fin :
		for line in fin :
			mxSize += 1
	#end with

	return mxSize
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
# Function: Read in a file as a line-by-line list of items
# Input ----
#   fname, str: path + name of the the sample files
# Returns ----
#   fItems, str list: ordered list of items from file
def readFileAsList(fname) :

	# ERROR CHECK: verify file exists
	verifyFile('', fname, False)

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
# Given a list of items, remove any items not in specified dict
# Input ----
#	theList, xxx list: list of items that may not be in theDict
#	theDict, xxx dict: dictionary against which to check (the keys)
# Returns ----
#   inList & outList, xxx list: which items from the
#       given list are found in dict & which aren't
def checkListAgainstDictKeys(theList, theDict) :

	# The items to return
	inList = list()
	outList = list()

	# extract the keys as a set
	keySet = set(theDict.keys())

	# Sift through the sample
	for item in theList :
		if item in keySet :
			inList.append(item)
		else :
			outList.append(item)
		#end if
	#end loop

	inList.sort()
	outList.sort()
	return inList, outList
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
# Function: save the given matrix as a .npy file
# Input ----
#	matrix, (NxN) list: the values to save
#	mname, str: name of the file to save
#	mpath, str: path to the folder to save the file
#	integer, bool: True means save values as int()
# Returns ----
#	nothing
def saveMatrixNumpy(matrix, mname, mpath, integer) :

	# Whether to save matrices as .txt or compressed .gz
	matrixExt = '.gz'

	# If folder doesn't exist, create it
	if not os.path.exists(mpath) :
		os.makedirs(mpath)
	#end if

	# Write to the file
#TODO: check if the fmt option is appropriate, or pointless
#	np.save(mpath + mname, matrix)
	if integer :
		np.savetxt(mpath + mname + matrixExt, matrix, fmt='%u')
	else :
		np.savetxt(mpath + mname + matrixExt, matrix, fmt='%f')
	#end if

#NOTE: In this case, the text file from savetxt() is much
#	smaller than the binary file from save()

#	# VERIFICATION: save as a text-readable file
#	if saveTextCopy :
#		saveMatrixText(matrix, "t"+mname, mpath, integer)
#	#end if

	return
#end def ######## ######## ######## 