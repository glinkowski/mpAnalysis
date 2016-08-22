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
# # select only meta-paths of specific length(s)
# limitMPLen = [1, 2, 3]
# 	# an empty list results in using all mp
# File name containing path z-score vectors
fZScoreSim = 'features_ZScoreSim.gz'

# verbose feedback ?
verbose = True

#TODO: remove / hardcode the following parameters


# Adjustible classifier parameters
useCfier = 1
	# 1 = Lasso, 2 = ElasticNet, ...
usePos = True
	# True/False: limit to only Positive coefficient values
useFeatPathSim = False
	# True/False: use the pathsim sum features
useFeatPathZScore = True
	# True/False: use the pathsim sum features
useFeatTermWeights = False
	# True/False: use the indirect term features
useFeatNeighbor = False
	# True/False: use the neighborhood features
useGivenRange = np.linspace(0.00005, 0.002, num=17)
	# array of vals; 'None' means to auto-search for alphas
maxClusters = 11
	# maximum number of clusters to use
	#	-1 = no maximum value is set


# LASSO params
lMaxIter = 500
lNorm = True
lFitIcpt = True


# text delimiter in output files
textDelim = '\t'

####### ####### ####### ####### 




####### ####### ####### ####### 
# ANCILLARY FUNCTIONS

def getGeneIndexLists(path, gDict) :

	# Create index lists for Known, Hidden, Unknown, TrueNeg
	gKnown = mp.readFileAsList(path + 'known.txt')
	giKnown = mp.convertToIndices(gKnown, geneDict)
	gHidden = mp.readFileAsList(path + 'concealed.txt')
	giHidden = mp.convertToIndices(gHidden, geneDict)
	giUnknown = [g for g in geneDict.values() if g not in giKnown]
	giTrueNeg = [g for g in giUnknown if g not in giHidden]

	return giKnown, giUnknown, giHidden, giTrueNeg
#end def ####### ####### ####### 


####### ####### ####### ####### 




####### ####### ####### ####### 
# PRIMARY FUNCTION

def predictIterative() :
	tstart = time.time()
	print("\nPerforming regression(s) on {}".format(sDir))

	mp.setParamVerbose(verbose)


	# 0) Create the useLabel variable
	# string: label for the output files
	# ie: ClusVote_c<Las/Enet/Log/SVM><P for Pos>_f<P for pathsim><Z for z-score>
	#	<T for term weights><N for neighborhood>
	if maxClusters > 0 :
		useLabel = 'Clus{}Vote_c'.format(maxClusters)
	else :
		useLabel = 'ClusVote_c'
	#end if
	if useCfier == 1 :
		useLabel = useLabel + 'Las'
	elif useCfier == 2 :
		useLabel = useLabel + 'Enet'
	else :
		print("ERROR: useCfier value is unrecognized: {}".format(useCfier))
	#end if
	if usePos :
		useLabel = useLabel + 'Pos'
	useLabel = useLabel + '_f'
	if useFeatPathSim :
		useLabel = useLabel + 'P'
		if limitMPLen :
			for item in limitMPLen :
				useLabel = useLabel + '{}'.format(item)
	if useFeatPathZScore :
		useLabel = useLabel + 'Z'
	if useFeatTermWeights :
		useLabel = useLabel + 'T'
	if useFeatNeighbor :
		useLabel = useLabel + 'N'
	#end if


	# 1) Load the gene-index dictionary & path names
	geneDict, pathDict = mp.getGeneAndPathDict(sDir)
	geneNames = list(geneDict.keys())
	geneNames.sort()
	pathNames = mp.removeInvertedPaths(pathDict)
	del pathDict

	# limit the metapaths by length
	#	part 1: get the desired indices
	idx = -1
	featPSIdx = list()
	for name in pathNames :
		idx += 1
		pLen = name.count('-') + 1
		if pLen in limitMPLen :
			featPSIdx.append( int(idx) )
	#end loop
	print("Limiting metapaths to {}, of length(s) {}".format(
		len(featPSIdx), limitMPLen))


	# 2) Load the network general features
	numFN = 0
	if useFeatNeighbor :
		featNbVals, featNbNames = mp.getFeaturesNeighborhood(sDir, 'LogScale')
		featNbNames = np.ravel(featNbNames)
		numFN = len(featNbNames)
	#end if
	numTW = 0
	if useFeatTermWeights :
		featTWVals, featTWNames = mp.getFeaturesTerms(sDir, 'Orig')
		featTWNames = np.ravel(featTWNames)
		numTW = len(featTWNames)
	#end if


	# 3) Loop over the list of the sample subdirectories
	dSubDirs = mp.getSubDirectoryList(dRoot+dDir)

	thisRound = 0
	#for si in dSubDirs[0:1] :
	for si in dSubDirs :
		# Display directory to examine
		sv = si.split('/')
		print("\n{}/{}/".format(sv[-3],sv[-2]))

		# Create index lists for Known, Hidden, Unknown, TrueNeg from files
		giKnown, giUnknown, giHidden, giTrueNeg = getGeneIndexLists(si, geneDict)


		# 4) Load the sample-specific features
		# PathSim features
		numFP = 0
		if useFeatPathSim :
			featPSVals = np.loadtxt(si + fSimilarity)
			# NOTE: previous version of mpPredict04a had extra zeros at end of vectors;
			#	discarding those columns
			featPSVals = featPSVals[:,0:len(pathNames)]
			featPSNames = pathNames
			numFP = len(featPSNames)

			# limit the metapaths by length
			#	part 2: keep only the desired columns
			if limitMPLen :
				featPSVals = featPSVals[:,featPSIdx]
				newFeatPSNames = list()
				for idx in featPSIdx :
					newFeatPSNames.append(featPSNames[idx])
				featPSNames = newFeatPSNames
				numFP = len(featPSNames)
		#end if

		# z-score of path counts features
		if useFeatPathZScore :
			featZSVals = np.loadtxt(si + fZScoreSim)
			featZSVals = featZSVals[:,0:len(pathNames)]
			featZSNames = pathNames
			numFP = len(featZSNames)
		#end if


		# 5) Combine the features as specified by parameters (useFeat...)
		features = np.zeros( (len(geneDict), 0), dtype=np.float32)
		featNames = list()

		if useFeatPathSim :
			print("    ... including PathSim sum features")
			features = np.hstack( (features, featPSVals) )
			featNames.extend(featPSNames)
		if useFeatPathZScore :
			print("    ... including path z-score features")
			features = np.hstack( (features, featZSVals) )
			featNames.extend(featZSNames)
		if useFeatNeighbor :
			print("    ... including neighborhood features")
			features = np.hstack( (features, featNbVals) )
			featNames.extend(np.ravel(featNbNames))
		if useFeatTermWeights :
			# Remove terms with no connection to gene set
			sumFTV = np.sum(featTWVals[giKnown,:], axis=0)
			keepIdx = np.nonzero(sumFTV)
			numTW = len(keepIdx[0])
			print("    ... including term membership features")
			features = np.hstack( (features, featTWVals[:,keepIdx[0]]) )
			featNames.extend(np.ravel(featTWNames[keepIdx]))
		# verify some features have been loaded
		if features.shape[1] == 0 :
			print("ERROR: No features were specified for classification.")
			sys.exit 
		#end if

		# Normalize the feature values
		features = mp.normalizeFeatureColumns(features)




		# 6) Prepare the test/train vectors & labels
		
#		# Create the structure to rank the Unknown genes & paths
#		gRankSums = np.zeros( len(geneList) )
#		pValSums = np.zeros( len(pathNames) )





#end def ####### ####### ####### 




####### ####### ####### ####### 
# FUNCTION CALL