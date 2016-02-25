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
# Function: read in the first column of a file
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
	count = 1
	firstWarn = True
	lvLen = -1
	for line in theFile :
		# Skip first nSkip number of rows
		if count <= nSkip :
			continue
		#end if

		line = line.rstrip()
		lv = line.split(textDelim)

		# ERROR CHECK: verify iCol is within range
		if (firstWarn == True) and (iCol >= len(lv)) :
			print( "WARNING: File contains {} columns, but col " +
				" {} was requested.".format(len(lv), iCol) )
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