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



####### ####### ####### ####### 
# PARAMETERS

# Variables/Quantities
#percHide = [0, 10, 25, 33, 50]		# percent of genes to conceal


# Input names & locations
useNtwk = 1		# network & samples to use (0 means fake)
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


# verbose feedback ?
newVerbose = True

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()
print("")


mp.setParamVerbose(newVerbose)


# load the genedict
# get the sample/folder list
# for each sample,
#   get, convert known genes to indices
#   perform lasso on path stats
#   output top paths & weights

# rank genes

# compare the two metrics

