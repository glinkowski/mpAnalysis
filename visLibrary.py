# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# functions used in creating visualizations
#
# These functions were created to aid in the creation of
#	visualizations used while analyzing the data and
#	results.
#
# Functions provided:
#	setParamVerbose(newVal)
#	setParamTextDelim(newVal)
#	setParamFileZeroPad(newMval, newOval)
#	readFileColumnAsString(fname, iCol, nSkip)
#	randSelectWithExclude(itemList, excludeSet, quantity)
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


# Length to pad the matrix file names:
fnMatrixZPad = 5
# Length to pad the output file names:
fnOutputZPad = 3
# Data delimiter to use in the output file:
textDelim = '\t'
# Whether to print non-error messages within these funcs
verbose = True

######## ######## ######## ######## 




######## ######## ######## ########
# Functions to set the global library parameters
# Input ----
#	NOTE: an input value of...
#		-1 keeps the parameter unchanged
#		-2 resets parameter to default
def setParamVerbose(newVal) :

	if newVal :
		verbose = True
	else :
		verbose = False
	#end if

	return
#end def ######## ######## ######## 
def setParamTextDelim(newVal) :

	if str(newVal == '-1') :
		textDelim = textDelim
	elif str(newVal == '-2') :
		textDelim = '\t'
	else :
		textDelim = str(newVal)
	#end if

	return
#end def ######## ######## ######## 
def setParamFileZeroPad(newMval, newOval) :

	if str(newMval == -1) :
		fnMatrixZPad = fnMatrixZPad
	elif str(newMval == -2) :
		fnMatrixZPad = 5
	else :
		fnMatrixZPad = newMval
	#end if
	if str(newOval == -1) :
		fnOutputZPad = fnOutputZPad
	elif str(newOval == -2) :
		fnOutputZPad = 3
	else :
		fnOutputZPad = newMval
	#end if

	return
#end def ######## ######## ######## 




######## ######## ######## ########
# Function: read in the desired column of a file
#	(typically a list of genes), skip N header rows
# Input ----
#	fname, str: path & filename of file to read
#	iCol, int: index of column to read in (! INDEXED AT 0 !)
#	nSkip, int; number of rows to skip at the top
# Returns ----
#	theList, str list: list of the items in the 
def readFileColumnAsString(fname, iCol, nSkip) :

	# ERROR CHECK: verify file exists
	if not os.path.isfile(fname) :
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	# the item to return
	theList = list()

	# Read in from the file
	theFile = open(fname, 'rb')
	count = 0
	firstWarn = True
	lvLen = -1
	for line in theFile :
		# Skip first nSkip number of rows
		if count < nSkip :
			continue
		#end if

		line = line.rstrip()
		lv = line.split(textDelim)

		# ERROR CHECK: verify iCol is within range
		if (firstWarn == True) and (iCol >= len(lv)) :
			print( "WARNING: File contains {} columns,".format(len(lv)) +
				" but col {} was requested.".format(iCol) )
			print( "    Note: Indexing starts at 0.")
			print( "    Returning last column. ")
			iCol = len(lv) - 1
			firstWarn = False
		#end if

		# ERROR CHECK: Warn if column length changes
		if (lvLen != len(lv)) and (lvLen >= 0) :
			print( "WARNING: Inconsistent number of columns in file." )
			print( "    Previous width: {}; Current: {}".format(
				lvLen, len(lv)) )
		#end if
		lvLen = len(lv)

		theList.append( str(lv[iCol]) )
		count += 1
	#end loop
	theFile.close()


	return theList
#end def ######## ######## ######## 




######## ######## ######## ########
# Function: randomly select items from a list, ensuring
#	none of the chosen items are in the exclude set
# Input ----
#	itemList, XX list: list of items to select from
#	excludeSet, XX set: set of items to exclude
#	quantity, int: number of items to select
# Returns ----
#	theList, XX list: list of randomly chosen items
def randSelectWithExclude(itemList, excludeSet, quantity) :

	# ERROR CHECK: if (quantity + # excludes) > # items,
	#	then return all items, less the excludes
	if quantity > ( len(itemList) - len(excludeSet) ) :
		print ( "WARNING: Improperly defined quantity ..." )
		print ( "    main item list contains {}".format(len(itemList)) )
		print ( "    exclude set contains {}".format(len(excludeSet)) )
		print ( "    desired size of selection {}".format(quantity) )
		theList = [x for x in itemList[:] if x not in excludeSet]
		print ( "    Returning {} items.".format(len(theList)) )
		return theList
	#end if


	# the item to return
	theList = [x for x in itemList[:] if x not in excludeSet]
	theList = random.sample(theList, quantity)
	theList.sort()

	return theList
#end def ######## ######## ######## 




######## ######## ######## ########
# Function: read in the top N items from a file
# Input ----
#	fname, str: path & filename of file to read
#		FORMAT: Col 0 is score/rank, Col 1 is the items
#	nItems, int: number of paths to return
#	nSkip, int; number of rows to skip at the top
# Returns ----
#	topNPaths, str list: list of top N paths, ordered by rank
#	rankScore, int? list: list of corresponding scores or ranks
def getTopRankedItems(fname, nItems, nSkip) :

	# ERROR CHECK: verify file exists
	if not os.path.isfile(fname) :
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	# The items to return
	topNPaths = list()
	rankScore = list()


	# Read in from the file
	theFile = open(fname, 'rb')
	count = 0
	for line in theFile :
		# Skip first nSkip number of rows
		if count < nSkip :
			continue
		# Exit if all desired items have been collected
		if count >= nItems :
			break
		#end if

		line = line.rstrip()
		lv = line.split(textDelim)

		# ERROR CHECK: verify expected columns exist
		if len(lv) < 2 :
			print ( "ERROR: The file is expected to contain at least " +
				"2 columns; {}".format(fname) )
			sys.exit()
		#end if

		topNPaths.append(lv[1])
		rankScore.append(lv[0])
		count += 1
	#end loop
	theFile.close()


	# ERROR CHECK: Warn if less items exist than asked for
	if count < nItems :
		print ( "WARNING: {} items were requested, but only ".format(count) +
			"{} exist in {}.".format(nItems, fname) )
	#end if

	return topNPaths, rankScore
#end def ######## ######## ######## 




######## ######## ######## ########
# Function: count the number of lines in a file
# Input ----
#	fname, str: path & filename of file to read
# Returns ----
#	cRows, int: total # of lines in the file
#	cColMin, int: minimum # of columns in the file
#	cColMax, int: maximum # of columns in the file
def countLinesInFile(fname) :

	# ERROR CHECK: verify file exists
	if not os.path.isfile(fname) :
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	cRows = 0
	cColMin = -1
	cColMax = 0

	if fname[-3:0] == '.gz' : 
		with gzip.open(fname, 'rb') as fin :
			for line in fin :
				lv = line.split(textDelim)
				lineLength = len(lv)

				cRows += 1
				
				if cColMax < lineLength :
					cColMax = lineLength

				if cColMin == -1 :
					cColMin = lineLength
				elif cColMin > lineLength :
					cColMin = lineLength
			#end loop
		#end with
	else :
		fin = open(fname, 'rb')
		for line in fin :
			lv = line.split(textDelim)
			lineLength = len(lv)

			cRows += 1
			
			if cColMax < lineLength :
				cColMax = lineLength

			if cColMin == -1 :
				cColMin = lineLength
			elif cColMin > lineLength :
				cColMin = lineLength
		#end loop
		fin.close()
	#end if

	return cRows, cColMin, cColMax
#end def ######## ######## ######## 




######## ######## ######## ########
# Function: calculate Recall (TPR), FPR, & Precision
# Input ----
#	path, str: path to file to read
#	name, str: name of file to read
# Returns ----
#	recall, 
#	FPR, 
#	precision, 
#	nHidden
def getAUCstats(path, name) :

	# Read in the ranked genes
	gHidden = readFileColumnAsString(path+'concealed.txt', 0, 0)
	gHidSet = set(gHidden)
	nHidden = len(gHidden)

	# In case concealed.txt is empty
	if nHidden == 0 :
		print("There are no Concealed Positives to predict in {}".format(path))
#		noHid = list()
#		return noHid, noHid, noHid, noHid
		return [-1], [-1], [-1], [-1]
	#end if

	# Declare the confusion matrix
	rows, colMin, colMax = countLinesInFile(path+name)
#	rows, colMin, colMax = countLinesInFile(path+'ranked_genes.txt')
	posActual = len(gHidden)
	negActual = rows - posActual
	confusion = np.zeros([2,rows])	# TP, FP


	# Create the matrix by reading in the file
	#	matrix is running count of: True Positives, False Positives
	fin = open(path+name)
#	fin = open(path+'ranked_genes.txt')
	col = 0
	TP = 0
	FP = 0
	for line in fin :
		line = line.rstrip()
		lv = line.split(textDelim)

		if lv[1] in gHidSet :
			TP += 1
		else :
			FP += 1
		#end if
		confusion[:,col] = [TP, FP]

		col += 1
	#end loop
	fin.close()


	# Convert into FPR, Recall (TPR), & Precision
	# TPR = TP / allP = confusion[0,i] / posActual
	# FPR = FP / allN = confusion[1,i] / negActual
	# recall = TPR
	# precision = TP / (TP + FP) = confusion[0,i] / (confusion[0,i] + confusion[1,i])
	recall = confusion[0,:] / posActual		# aka TPR
	FPR = confusion[1,:] / negActual
	precision = confusion[0,:] / (confusion[0,:] + confusion[1,:])


	return FPR, recall, precision, nHidden
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: choose an unused name for the output file
# Input ----
#   path, str: path to the network files
#   name, str: name of the network to use
#	extension, str: extension to append to file
# Returns ----
#   fname, str: name of output file (without path)
def nameOutputFile(path, name, extension) :

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
	fname = name + "-{}.".format(str(num).zfill(zpad)) + extension
#	fname = name + "-{}.txt".format(str(num).zfill(zpad))
	while fname in fileSet :
		num += 1
		fname = name + "-{}.".format(str(num).zfill(zpad)) + extension
	#end loop

	return fname
#end def ######## ######## ######## 
