# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Approach 2.5: learn top paths from PathSim sum, find genes
#	Version 2, Step 2
#	iteratively predict set membership
#		Lasso Pos, 1-Class sampling, small vote
#		features include PathSim, neighborhood, terms
#	Requires: mpPredict04a
#
# This approach finds similarity of each gene to a set
#	along each metapath, then learns the most important
#	paths, which are treated as features.
# Then, the prediction of set membership is performed
#	iteratively, allowing the method to separate out
#	the more difficult items to predict.
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
sDir = '../Dropbox/mp/output/pred04-msig100'

# File name containing pathsim feature vectors
fSimilarity = 'features_PathSim.gz'
# select only meta-paths of specific length(s)
limitMPLen = [1, 2, 3]
	# an empty list results in using all mp

# verbose feedback ?
newVerbose = True

#TODO: remove / hardcode the following parameters


# Adjustible classifier parameters
useCfier = 1
	# 1 = Lasso, 2 = ElasticNet, ...
usePos = True
	# True/False: limit to only Positive coefficient values
useFeatPaths = True
	# True/False: use the pathsim sum features
useFeatNeighbor = False
	# True/False: use the neighborhood features
useGivenRange = np.linspace(0.00005, 0.002, num=17)
	# array of vals; 'None' means to auto-search for alphas
maxClusters = 11
	# maximum number of clusters to use
	#	-1 = no maximum value is set


# LASSO params
lMaxIter = 1000
lNorm = True
lFitIcpt = True


# text delimiter in output files
textDelim = '\t'

####### ####### ####### ####### 




####### ####### ####### ####### 
# ANCILLARY FUNCTIONS

####### ####### ####### ####### 




####### ####### ####### ####### 
# PRIMARY FUNCTION

def predictIterative() :
	tstart = time.time()
	print("\nPerforming regression(s) on {}".format(sDir))

	mp.setParamVerbose(newVerbose)
#end def ####### ####### ####### 




####### ####### ####### ####### 
# FUNCTION CALL