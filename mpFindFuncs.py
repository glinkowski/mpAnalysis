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
#	readGeneFile(fname)
# ---------------------------------------------------------

import os.path
import sys



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
def readGeneFile(fname) :

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