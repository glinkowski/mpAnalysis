# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# clean out files created by mpPredict04 process
#	not including features_PathSim.gz
#
# This cleans out files created when testing the code.
# The following will be deleted from each subdirectory:
#	AUC-(...).png
#	parameters-(...).txt
#	ranked_features_Top(N)-(...).txt
#	ranked_genes-(...).txt
#	ranked_paths-(...).txt
# ---------------------------------------------------------

import mpLibrary as mp
import os



####### ####### ####### ####### 
# PARAMETERS

# folder containing the samples & results
dDir = 'pred04-achilles200'
dRoot = '../Dropbox/mp/output/'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION


print("\nCleaning mpPredict04 files from {}".format(dDir))

# Get a list of the subdirectories
if not dDir.endswith('/') :
	dDir = dDir + '/'
if not dRoot.endswith('/') :
	dRoot = dRoot + '/'
dSubDirs = mp.getSubDirectoryList(dRoot + dDir)


# Delete specified files from each subdir
for si in dSubDirs :

	fileList = os.listdir(si)
	fileList.sort()

	for fn in fileList :
		if fn.startswith('parameters-') and fn.endswith('.txt') :
			os.remove(si + fn)
		elif fn.startswith('ranked_features_Top') and fn.endswith('.txt') :
			os.remove(si + fn)
		elif fn.startswith('ranked_genes-') and fn.endswith('.txt') :
			os.remove(si + fn)
		elif fn.startswith('ranked_paths-') and fn.endswith('.txt') :
			os.remove(si + fn)
#		elif fn.startswith('AUC-') and fn.endswith('.png') :
#			os.remove(si + fn)
	#end loop
#end loop



print("\nDone.\n")