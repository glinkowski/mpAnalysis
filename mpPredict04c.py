# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Approach 2: learn top paths from PathSim sum, find genes
#	Version 2, Step 2(c)
#	LASSO + Voting
#
# The first approach tries to first find the most
#	important metapaths, then rank genes by similarity
#	along those. This approach instead finds the set
#	similarity along each metapath, then tries to learn
#	the most important paths.
#
# Step 2c: From the pre-built per-gene feature vector,
#	apply LASSO linear regression to select a small number
#	of important paths. Follow these paths to find genes
#	connected to the set, and rank those by similarity.
# NEW: Do this several times, and vote among the gene
#	rankings to create a final rank
# ---------------------------------------------------------

import mpLibrary as mp
import time
import numpy as np
import gzip
from sklearn import linear_model as lm
import random



####### ####### ####### ####### 
# PARAMETERS

# folder containing the pre-processed samples
dDir = 'pred04-batch-000'

# Input names & locations
useNtwk = 0		# network & samples to use (0 means fake)
if useNtwk == 0 :
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
	dRoot = 'outputFake/'
else :
	eName = 'all_v1_g2e11t0'
	ePath = '../Dropbox/mp/networks/'
	dRoot = '../Dropbox/mp/output/'
#end if

# File name containing feature vectors
fSimilarity = 'features_PathSim.gz' 


# LASSO params
lAlpha = 0.01
lMaxIter = 10000
lNorm = True
lPos = False
lFitIcpt = True
lSelctn = 'random' # random vs cyclic

pLabel = 1
nLabel = 0

# verbose feedback ?
newVerbose = True


textDelim = '\t'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()
print("")

mp.setParamVerbose(newVerbose)



# 1) Load the gene-index dictionary & path names
print("Creating the gene-index dictionary.")
geneDict = mp.readGenesFile(ePath, eName)
geneList = list(geneDict.keys())
geneList.sort()

print("Reading in the path names.")
pathDict = mp.readKeyFile(ePath, eName)
pathNames = mp.removeInvertedPaths(pathDict)
del pathDict



# 2) Get a list of the sample subdirectories
dSubDirs = mp.getSubDirectoryList(dRoot+dDir)
