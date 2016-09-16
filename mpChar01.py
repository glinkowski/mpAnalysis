# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Characterization: define unique types of connections in set
#	Step 1:
#		create z-score features
#
# For each sample, perform once on the whole set and then
#	once on each of 4 folds. Partition set into Known Pos
#	and Hidden Pos.
# Then, create the metapath z-score feature, a per-gene
#	feature vector where each entry is the sum of that
#	gene's similarity to each member of the Known Pos set.
# ---------------------------------------------------------

import charLibrary as cl
import time
import numpy as np
import gzip



####### ####### ####### ####### 
# PARAMETERS

# Input names & locations
useNtwk = 0		# network & samples to use (0 means fake)
if useNtwk == 0 :
#	eName = 'fakeNtwk00_g2e3t10'
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
	sPath = 'samplesFake/'
	oRoot = 'outputFake/'
else :
	eName = 'all_v3beta_g2e9t0'
	ePath = '../Dropbox/mp/networks/'
	sPath = '../Dropbox/mp/samples/MSIGDB/select01/'
#	sPath = '../Dropbox/mp/samples-test1/'
	oRoot = '../Dropbox/mp/output/'
#end if


# verbose feedback ?
verboseOutput = False

####### ####### ####### ####### 




####### ####### ####### ####### 
# ANCILLARY FUNCTIONS

#end def ####### ####### ####### 

####### ####### ####### ####### 




####### ####### ####### ####### 
# PRIMARY FUNCTION

def createFeatureZScore(eDir, sDir, oDir, printFlag) :
# Partition the samples and create z-score features
# Input ----
#	eDir, str: directory containing edge list
#	sDir, str: directory containing samples
#	oDor, str: output directory
#	printFlag, bool: whether to enable terminal output
# Returns ----
#	nothing
# Creates ----
#	directory in oDir containing folder for each sample
#	(whole + 4 folds), containing gene partitions, and
#	a matrix of z-score feature vectors for that sample



#end def ####### ####### ####### 




####### ####### ####### ####### 
# FUNCTION CALL

tstart = time.time()
# print("\nPredicting gene ranks for {}".format(sDir))

# # Main Function call
# print("Calling function predictIterative() ...")
# predictIterative(True)

# print("--elapsed time: {:.3} (s)".format(time.time()-tstart))



print("\nDone.\n")