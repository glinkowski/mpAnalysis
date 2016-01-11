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
# ---------------------------------------------------------

import os.path
import sys
import numpy as np




######## ######## ######## ######## 
# Function: Read in a file as a line-by-line list of items
# Input:
#   fname, str - path + name of the the sample files
# Returns:
#	fItems, list of str - ordered list of items from file
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
# Input:
#	fname, str - path & name to keep file
# Returns:
#	gDict, dict - 
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
	gDict = dict()
	gf = open(fname, "rb")
	index = 0
	for line in gf :
		gene = line.rstrip()	# remove "\n"
		gDict[gene] = index
		index += 1
	#end loop

	return gDict
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Read in the dataset from a samplename
#	Check for variants: ".txt", "_UP.txt", "_DN.txt"
# Input:
#   fname, str - path + name of the the sample files
#   up, down: boolean -- only read the _UP or _DN file if true
# Returns:
#	sNodes, list of str - ordered list of names from file(s)
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
# Input:
#	path, str - path to the network files
#	name, str - name of the network to use
# Returns:
#	mpDict, dict
#		key, str: name of metapath
#		value, tuple: matrix/file ID number
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
	fn = open(fname, "rb")
	firstline = True
	for line in fn :

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

	return keyDict
#end def ######## ######## ######## 



