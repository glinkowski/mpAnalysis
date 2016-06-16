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
import time



####### ####### ####### ####### 
# PARAMETERS

# The network to use and directory path
#eName = 'all_v3beta_g2e9t0'
#ePath = '../Dropbox/mp/networks/'
eName = 'fakeNtwk00_g2e3t10'
ePath = 'networks/'

####### ####### ####### ####### 


def reverseName(pName) :
	pv = pName.split('-')
	rv = pv[::-1]
	revName = ''
	for part in rv :
		revName = revName + part + '-'
	revName.rstrip('-')

	return revName
#end def ############### 


####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()

print("\nCreating neighborhood features for {}".format(eName))


# 1) Load the primary matrices to memory
print("Loading the primary path matrices ...")
if not ePath.endswith('/')
	ePath = ePath + '/'
if not eName.endswith('/')
	eName = eName + '/'
pathDict = pp.readKeyFilePP(ePath + eName)
pathNames = list(pathDict.keys())
pathNames.sort()

print(pathNames)



# 2) Allocate space for the features
nCols = 0
for name in pathNames :
	nCols += 1
	revName = reverseName(name)
	if revName != name :
		nCols += 1
#end loop

nRows = -1
fName = (ePath + eName.rstrip('/') + "_MetaPaths/" +
			"{}.gz".format(str(0).zfill(fnMatrixZPad)) )
with gzip.open(fname, 'rb') as fin :
	for line in fin :
		nRows += 1
#end with

featVals = np.zeros( (nRows,nCols), dtype=np.float32 )
featNames = np.zeros( nCols, dtype=object )
