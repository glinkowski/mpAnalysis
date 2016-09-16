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
import random
import gzip
from sklearn import cluster as skc
import warnings




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
#end def ####### ####### ####### 