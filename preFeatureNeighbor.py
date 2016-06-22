# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# Pre-Processing of the network
#	Create feature vectors of each gene's neighborhood
#
# This assumes the pre-processing has already been run.
#	From the path-count matrices, create a feature vector
#	describing the gene's neighborhood.
# ---------------------------------------------------------

import preProcFuncs as pp
import matplotlib.pyplot as plt
import numpy as np
import gzip
import time
import mpLibrary as mp



####### ####### ####### ####### 
# PARAMETERS

# The network to use and directory path
eName = 'all_v3beta_g2e9t0'
ePath = '../Dropbox/mp/networks/'
#eName = 'fakeNtwk00_g2e3t10'
#ePath = 'networks/'

# expected size of matrix file name
keyZPad = 5
####### ####### ####### ####### 


def reverseName(pName) :
	pv = pName.split('-')
	rv = pv[::-1]
	revName = ''
	for part in rv :
		revName = revName + part + '-'
#	revName.rstrip('-')
	revName = revName[0:-1]

	return revName
#end def ############### 
def concatenatePaths(path01, path02) :
	if not path01.endswith('/') :
		path01 = path01 + '/'
	if not path02.endswith('/') :
		path02 = path02 + '/'

	return (path01 + path02)
#end def ############### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()

print("\nCreating neighborhood features for {}".format(eName))


# 1) Load the primary matrices to memory
print("Loading the primary path matrices ...")
fPath = concatenatePaths(ePath, eName)
pathDict = pp.readKeyFilePP(fPath)
#pathNames = list(pathDict.keys())
pathNames = mp.removeInvertedPaths(pathDict)
pathNames.sort()

#print(pathNames)
print("len pathNames = {}; unique path names = {}".format( len(pathNames), len(np.unique(pathNames)) ))



# 2) Allocate space for the features

# Determine number of columns
nCols = 0
for name in pathNames :
	nCols += 1
	revName = reverseName(name)
	if revName != name :
		nCols += 1
#end loop
# nCols = len(pathNames)

# Determine number of rows
nRows = 0
mpPath = ePath + eName.rstrip('/') + "_MetaPaths/"
fName = (mpPath + "{}.gz".format(str(0).zfill(keyZPad)) )
with gzip.open(fName, 'rb') as fin :
	for line in fin :
		nRows += 1
#end with

featVals = np.zeros( (nRows,nCols), dtype=np.float32 )
featNames = np.zeros( nCols, dtype=object )



# 3) Load each path-count matrix 
print("Creating features from the metapath count matrices ...")
pathList = list()
col = -1
mcount = 0
for name in pathNames :

	mcount += 1
	if not (mcount % 25) :
		print("  reading matrix # {}: {}".format( mcount, name ))

	if name not in pathList :
		col += 1
#		if not (col % 25) :
#			print("  saving feature # {}: {}".format(col, name))

		matrix = pp.getPathMatrix(pathDict[name], mpPath, '', nRows)

		# create the neighborhood feature
		featVals[:,col] = np.sum(matrix, axis=1) - np.diag(matrix)
		pathList.append(name)

		# do the same for the inverse path
		nameRev = reverseName(name)
		if nameRev not in pathList :
			col += 1
#			if not (col % 25) :
#				print("  saving feature # {}: {}".format(col, name))
			featVals[:,col] = np.sum(matrix, axis=0) - np.diag(matrix)
			pathList.append(nameRev)
		#end if
#end loop



# 4) Apply log() to fix scaling of feature values
print("Applying logarithmic scaling to the feature values...")
featMod = np.add(featVals, 1.0)
featMod = np.log(featMod)



# 5) Save the features to file 
print("Saving the feature matrices to file...")
fPath = concatenatePaths(ePath, eName)

# write the feature names
with open(fPath + 'featNeighbor_Names.txt', 'w') as fout :
	fout.write('NL{}-{}'.format( (pathList[0].count('-') + 1), pathList[0] ))
	for i in range(1, len(pathList)) :
		fout.write('\nNL{}-{}'.format( (pathList[i].count('-') + 1), pathList[i] ))
#end with

# write the feature matrices, original and modified
pp.setParamSaveTextCopy(True)
pp.saveMatrixNumpy(featVals, 'featNeighbor_Orig', fPath, False)
pp.saveMatrixNumpy(featMod, 'featNeighbor_LogScale', fPath, False)

print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))



print("\nDone.\n")