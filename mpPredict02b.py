# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Approach 1: path counting & percentile beat
#	Version 2, Step 2
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
import random


# For testing & verification
#import random
#random.seed(42)



####### ####### ####### ####### 
# PARAMETERS

# save location of the ranked sample data
sFolder = 'pred02-batch-012'
#sRoot = '../Dropbox/mp/output/'

# Number of top paths to use
numTopK = 10


# Input names & locations
useNtwk = 0		# network & samples to use (0 means fake)
if useNtwk == 0 :
#	eName = 'fakeNtwk00_g2e3t10'
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
	sRoot = '../Dropbox/mp/outputFake/'
	retCutoffs = [2, 3, 5, 9, 13]
else :
	eName = 'all_v1_g2e11t0'
	ePath = '../Dropbox/mp/networks/'
	sRoot = '../Dropbox/mp/output/'
	retCutoffs = [50, 100, 200, 500, 1000]
#end if


# verbose feedback ?
verbose = False

matrixDT = np.float32

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()
print("")

mp.setParamVerbose(verbose)
mp.setParamMatrixDT(matrixDT)



# 1) Load the gene-index dictionary & path names
print("Creating the gene-index dictionary.")
geneDict = mp.readGenesFile(ePath, eName)
geneList = list(geneDict.keys())
geneList.sort()

print("Reading in the path names.")
pathDict = mp.readKeyFile(ePath, eName)
pathNames = mp.removeInvertedPaths(pathDict)
#del pathDict



# 2) Get a list of the sample subdirectories
if not sRoot.endswith('/') :
	sRoot = sRoot + '/'
dSubDirs = mp.getSubDirectoryList(sRoot + sFolder)
#print(dSubDirs[0])



# 3) Read in the Ignore List (wich paths to ignore)
# TODO: this.
pathIgnore = []



# For each sample (each subdirectory), rank genes
#		save top genes & rank to a file
print("Ranking genes for each sample ...")
#for si in dSubDirs[1:2] :
for si in dSubDirs :

	# Display directory to examine
	sv = si.split('/')
	print("  {}/{}/".format(sv[-3],sv[-2]))

	# Read in the ranked paths
	pathRanked = mp.readRankedPaths(si)

	# Sort array DESCENDING by rank and ASCENDING by length
	pathRanked = np.sort(pathRanked, order=['stat_inverse', 'length', 'name'])

	# Read in the non-hidden genes
	# Create index lists for Known, Hidden, Unknown
	gKnown = mp.readFileAsList(si+'known.txt')
	giKnown = mp.convertToIndices(gKnown, geneDict)
	gHidden = mp.readFileAsList(si+'concealed.txt')
	giHidden = mp.convertToIndices(gHidden, geneDict)
	giUnknown = [i for i in range(len(geneDict)) if i not in giKnown]



	# 4) Select top paths for naive, guided01, random tests

	# Select top K paths, naive method
	topNaive = list()
	topNaiveScore = list()
	for i in range(numTopK) :
		thisPath = pathRanked['name'][i]
		# skip if contains item in Ignore list
		for igPath in pathIgnore:
			if igPath in thisPath :
				continue
		#end if
		topNaive.append(thisPath)
		topNaiveScore.append(pathRanked['stat'][i])
	#end if

	# Select top K paths, first guided method
	topGuided = list()
	topGuidedScore = list()
	added = 0
	skipPaths = list()
#	skipPathsLen = list()
	for i in range(numTopK) :
		thisPath = pathRanked['name'][i]
#TODO: check the reverse path?
		# thisPathList = thisPath.split('-')[::-1]
		# thisPathRev = ''
		# for i in range(len(thisPathList))[::-1] :
		# 	thisPathRev = thisPathRev + '-' + thisPathList[i]
		# thisP
		# skip if contains item in Ignore list
		for igPath in pathIgnore:
			if igPath in thisPath :
				continue
		# skip if contains a path already added
		for sp in skipPaths :
			if sp in thisPath :
#			if (sp in thisPath) and ((skipPathsLen + 1) == pathRanked['length'][i]) :
				continue
#			if thisPath in sp :
#				continue
#TODO: better comparison? percent of path is similar?
		#end if
		topGuided.append(thisPath)
		topGuidedScore.append(pathRanked['stat'][i])
		added += 1
		# add to the skip list
		if pathRanked['length'][i] >= 2 :
			skipPaths.append(thisPath)
#			skipPathsLen.append(pathRanked['length'][i])
	#end loop

	# Select K paths at random (for comparison)
	topRandom = list()
	topRandomScore = list()
	randIdx = random.sample(range(len(pathRanked)), numTopK)
	for i in randIdx :
		topRandom.append(pathRanked['name'][i])
		topRandomScore.append(pathRanked['stat'][i])
	#end loop

	mp.writeChosenPaths(si, 'naive', topNaive, topNaiveScore)
	mp.writeChosenPaths(si, 'guided', topGuided, topGuidedScore)
	mp.writeChosenPaths(si, 'random', topRandom, topRandomScore)



	# 5) Rank genes by similarity along selected metapaths
	#	For each path, load pathsim matrix, sum similarity

	# For each path get gene similarity, naive
	simArrayNaive = np.empty([len(giUnknown), len(topNaive)], dtype=matrixDT)
	idx = 0
	for p in topNaive :
		simMatrix = mp.getSimMatrix( pathDict[p], ePath,
			eName, (len(giKnown) + len(giUnknown)) )
		simCols = simMatrix[:,giKnown]
		simArrayNaive[:,idx] = np.sum(simCols, axis=1)[giUnknown]
		idx += 1
	#end loop

	# second Naive, using scores as weights
	simArrayNaive02 = np.copy(simArrayNaive)
	for c in range(simArrayNaive.shape[1]) :
		simArrayNaive02[:,c] = np.multiply(simArrayNaive02[:,c], topNaiveScore[c])
	#end loop

	# For each path get gene similarity, first guided
	simArrayGuided = np.empty([len(giUnknown), len(topGuided)], dtype=matrixDT)
	idx = 0
	for p in topGuided :
		simMatrix = mp.getSimMatrix( pathDict[p], ePath,
			eName, (len(giKnown) + len(giUnknown)) )
		simCols = simMatrix[:,giKnown]
		simArrayGuided[:,idx] = np.sum(simCols, axis=1)[giUnknown]
		idx += 1
	#end loop

	# second Guided, using scores as weights
	simArrayGuided02 = np.copy(simArrayGuided)
	for c in range(simArrayGuided.shape[1]) :
		simArrayGuided02[:,c] = np.multiply(simArrayGuided02[:,c], topGuidedScore[c])
	#end loop

	# For each path get gene similarity, first guided
	simArrayRandom = np.empty([len(giUnknown), len(topRandom)], dtype=matrixDT)
	idx = 0
	for p in topRandom :
		simMatrix = mp.getSimMatrix( pathDict[p], ePath,
			eName, (len(giKnown) + len(giUnknown)) )
		simCols = simMatrix[:,giKnown]
		simArrayRandom[:,idx] = np.sum(simCols, axis=1)[giUnknown]
		idx += 1
	#end loop

	# second Random, using scores as weights
	simArrayRandom02 = np.copy(simArrayRandom)
	for c in range(simArrayRandom.shape[1]) :
		simArrayRandom02[:,c] = np.multiply(simArrayRandom02[:,c], topRandomScore[c])
	#end loop



	# 6) Write the ranked_genes files + chosen paths

	# gRanks to be used for voting method
	gRanks = np.empty([len(giUnknown), 6], dtype=object)
	gRanks[:,0] = mp.writeRankedGenes02(si, 'naive', simArrayNaive,
		geneDict, giKnown, gHidden, retCutoffs)
	gRanks[:,1] = mp.writeRankedGenes02(si, 'naive02', simArrayNaive02,
		geneDict, giKnown, gHidden, retCutoffs)
	gRanks[:,2] = mp.writeRankedGenes02(si, 'guided', simArrayGuided,
		geneDict, giKnown, gHidden, retCutoffs)
	gRanks[:,3] = mp.writeRankedGenes02(si, 'guided02', simArrayGuided02,
		geneDict, giKnown, gHidden, retCutoffs)
	gRanks[:,4] = mp.writeRankedGenes02(si, 'random', simArrayRandom,
		geneDict, giKnown, gHidden, retCutoffs)
	gRanks[:,5] = mp.writeRankedGenes02(si, 'random02', simArrayRandom02,
		geneDict, giKnown, gHidden, retCutoffs)



	# 7) Create a ranking based on voting

	# get ranks for each gene
	#	invert the value to be passed to the function
	rankVote = np.zeros( (len(geneDict), gRanks.shape[1]))
	for rank in range(gRanks.shape[0]) :
		for col in range(gRanks.shape[1]) :
			row = geneDict[gRanks[rank,col]]
			rankVote[row,col] -= rank
	#end loop
	rankVote = rankVote[giUnknown,:]

	# write to file
	mp.writeRankedGenes02(si, 'voting', rankVote[:,[0,2,4]],
		geneDict, giKnown, gHidden, retCutoffs)
	mp.writeRankedGenes02(si, 'voting02', rankVote[:,[1,3,5]],
		geneDict, giKnown, gHidden, retCutoffs)
	mp.writeRankedGenes02(si, 'votingAll', rankVote,
		geneDict, giKnown, gHidden, retCutoffs)

#end loop



print("\nDone.\n")