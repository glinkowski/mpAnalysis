# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# A two-step approach to set membership prediction, using
#	a batch approach. This file contains Step Two:
#	Rank the genes for the top selected meta-paths.
#
# After meta-paths have been ranked for a batch of samples,
#	select the top K paths from the ranked list and 
#	calculate gene similarity along those paths.
# ---------------------------------------------------------

import mpLibrary as mp
import time
import numpy as np
import gzip


# For testing & verification
#import random
#random.seed(42)



####### ####### ####### ####### 
# PARAMETERS

# save location of the ranked sample data
sFolder = 'pred02-batch-000'
sRoot = '../Dropbox/mp/output/'


# Input names & locations
useNtwk = 0		# network & samples to use (0 means fake)
if useNtwk == 0 :
#	eName = 'fakeNtwk00_g2e3t10'
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
else :
	eName = 'all_v1_g2e11t0'
	ePath = '../Dropbox/mp/networks/'
#end if


# verbose feedback ?
verbose = False

####### ####### ####### ####### 

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()
print("")

mp.setParamVerbose(verbose)


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
if not sRoot.endswith('/') :
	sRoot = sRoot + '/'
dSubDirs = mp.getSubDirectoryList(sRoot + sFolder)
print(dSubDirs[0])


# 3) Read in the Ignore List (wich paths to ignore)
# TODO: this.
pathIgnore = []


# 4) For each sample (subdir), rank genes
#		save top genes & rank to a file

for si in dSubDirs[0:1] :
#for si in dSubDirs :

	# Display directory to examine
	sv = si.split('/')
	print(si)
	print(sv)
	print("\n{}/{}/".format(sv[-3],sv[-2]))


	# Read in the ranked paths
	pathRanked = mp.readRankedPaths(si)

	print(pathRanked[0])











print("\nDone.\n")