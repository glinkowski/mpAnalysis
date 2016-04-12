# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
#
# A different approach to ranking meta-paths.
# In this case, the collected counts between genes and the
#	training set will be normalized and used to weight
#	the metapaths by predictability, or importance. This
#	can be throught of as using feedback from the PathSim
#	metric to rank genes in order to influence and automate
#	the decision of which paths to use for prediction.
#
# This is the second part of the process: 
# ---------------------------------------------------------

import mpLibrary as mp
import time
import numpy as np
import gzip
from sklearn import linear_model as lm



####### ####### ####### ####### 
# PARAMETERS

# Variables/Quantities
#percHide = [0, 10, 25, 33, 50]		# percent of genes to conceal


# Input names & locations
useNtwk = 0		# network & samples to use (0 means fake)

if useNtwk == 0 :
#	eName = 'fakeNtwk00_g2e3t10'
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
#	sPath = 'samplesFake/'
	dRoot = '../Dropbox/mp/outputFake/'
else :
	eName = 'all_v1_g2e11t0'
	ePath = '../Dropbox/mp/networks/'
#	sPath = '../Dropbox/mp/samples-test1/'
	dRoot = '../Dropbox/mp/output/'
#end if

# Output path
dDir = 'pred03-batch-000'


# File names for similarity metrics
#fGroupNorm = 'Pxy.gz'
fGroupNorm = 'Pxy-mod-norm.gz'
fOrigSum = 'SxySum.gz'


# LASSO params
lAlpha = 0.0001
lMaxIter = 1000
lNorm = False
lCopy = False


# verbose feedback ?
newVerbose = True

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()
print("")


mp.setParamVerbose(newVerbose)


# 1) Load the gene-index dictionary
print("Creating the gene-index dictionary.")
geneDict = mp.readGenesFile(ePath, eName)


# 2) Get a list of the sample subdirectories
dSubDirs = mp.getSubDirectoryList(dRoot+dDir)

# 3) For each sample (subdir), perform LASSO
#		save top paths & weights to a file

for si in dSubDirs[1:2] :
#for si in dSubDirs :

	# Read in the known genes
	print(si)
	gKnown = mp.readFileAsList(si+'known.txt')
	print(gKnown)
	giKnown = mp.convertToIndices(gKnown, geneDict)
	print(giKnown)

	# Read in the (modified) group PathSim matrix
	simGroup = mp.readFileAsMatrix(si, fGroupNorm)
	simOrig = mp.readFileAsMatrix(si, fOrigSum)
	print(simGroup.shape, simOrig.shape)
#	print(simGroup)
	print(np.amax(simGroup,axis=1)[0:12])
#	print(simOrig)
	print(np.amax(simOrig,axis=1)[0:12])

	# Extract the rows on which to perform LASSO
	targetGroup = simGroup[giKnown,:]
	targetOrig = simOrig[giKnown,:]
	print targetGroup[0:3,0:10]

	# Perform LASSO
	cfier = lm.Lasso(alpha=lAlpha)
	cfier.fit(targetGroup, np.ones(targetGroup.shape[0]))
	print(cfier.coef_.shape, cfier.coef_)



# load the genedict
# get the sample/folder list
# for each sample,
#   get, convert known genes to indices
#   perform lasso on path stats
#   output top paths & weights

# rank genes

# compare the two metrics

