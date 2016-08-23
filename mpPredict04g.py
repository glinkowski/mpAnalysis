# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Approach 2: learn top paths from PathSim sum, find genes
#	Version 2, Step 2(g)
#	1-Class sampling, cluster Unknown, ensemble prediction,
#		multiple types of classifiers
#		include PathSim & Neighborhood features
#
# The first approach tries to first find the most
#	important metapaths, then rank genes by similarity
#	along those. This approach instead finds the set
#	similarity along each metapath, then tries to learn
#	the most important paths.
#
# Step 2g: This is an adaptation of the voting approach,
#	where the Unknowns are clustered, and then a predictor
#	is trained on each of those. The final results include
#	votes from each predictor.
#
# Updates: option to include the Neighborhood features
#	alongside the PathSim ones
# ---------------------------------------------------------

import mpLibrary as mp
import time
import numpy as np
import gzip
from sklearn import linear_model as lm
from sklearn import svm as svm
import random



####### ####### ####### ####### 
# PARAMETERS

# folder containing the pre-processed samples
sDir = '../Dropbox/mp/output/pred04-test01'

# File name containing pathsim feature vectors
fSimilarity = 'features_PathSim.gz'
# select only meta-paths of specific length(s)
limitMPLen = [1, 2, 3]
	# an empty list results in using all mp
# File name containing path z-score vectors
fZScoreSim = 'features_ZScoreSim.gz'


# verbose feedback ?
newVerbose = True


# Adjustible classifier parameters
useCfier = 5
	# 1 = Lasso, 2 = ElasticNet, 3 = SVM, ...
usePos = True
	# True/False: limit to only Positive coefficient values
useFeatPathSim = True
	# True/False: use the pathsim sum features
useFeatPathZScore = False
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

# Elastic Net params
enRatios = [0.2, 0.4, 0.6, 0.75, 0.85, 0.9, 0.95]
enNAlphas = 11
enMaxIter = 400
enFitIncept = True
enNorm = True
enCopy = True

# SVM params
svmMaxIter = 700
svmKernel = 'rbf' # linear, poly, rbf, sigmoid


# text delimiter in output files
textDelim = '\t'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()
print("\nPerforming regression(s) on {}".format(sDir))

mp.setParamVerbose(newVerbose)



# 0) Create the useLabel variable
# string: label for the output files
# ie: ClusVote_c<Las/Enet/Log><P for Pos>_f<P for pathsim><N for neighborhood>
if maxClusters > 0 :
	useLabel = 'Clus{}Vote'.format(maxClusters)
else :
	useLabel = 'ClusVote'
#end if
useLabel = useLabel + '_c'
if useCfier == 1 :
	useLabel = useLabel + 'Las'
elif useCfier == 2 :
	useLabel = useLabel + 'Enet'
elif useCfier == 3 :
	useLabel = useLabel + 'SVM'
elif useCfier == 4 :
	useLabel = useLabel + 'LinSVM'
elif useCfier == 5 :
	useLabel = useLabel + 'LasSVM'
else :
	print("ERROR: useCfier value is unrecognized: {}".format(useCfier))
#end if
if usePos and (useCfier < 3) :
	useLabel = useLabel + 'Pos'
#end if
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

print(useLabel)


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
# if newVerbose :
# 	print(featPSIdx)
# #end if



# 2) Load the network general features
numFN = 0
if useFeatNeighbor :
	featNbVals, featNbNames = mp.getFeaturesNeighborhood(sDir, 'LogScale')
	featNbNames = np.ravel(featNbNames)
	numFN = len(featNbNames)
#end if

numFT = 0
if useFeatTermWeights :
	featTWVals, featTWNames = mp.getFeaturesTerms(sDir, 'Orig')
	featTWNames = np.ravel(featTWNames)
	numFT = len(featTWNames)
#end if


# 3) Loop over all the subdirectories
dSubDirs = mp.getSubDirectoryList(sDir)

thisRound = 0
#for si in dSubDirs[0:1] :
for si in dSubDirs :

	# Display directory to examine
	sv = si.split('/')
	print("\n{}/{}/".format(sv[-3],sv[-2]))



	# 4) Load the PathSim features
	numFP = 0
	if useFeatPathSim :
#		featPSVals = mp.readFileAsMatrix(si, fSimilarity)
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

	# Load the ZScore features
	numFZ = 0
	if useFeatPathZScore :
#		featPSVals = mp.readFileAsMatrix(si, fSimilarity)
		featZSVals = np.loadtxt(si + fZScoreSim)
		# NOTE: previous version of mpPredict04a had extra zeros at end of vectors;
		#	discarding those columns
		featZSVals = featZSVals[:,0:len(pathNames)]
		featZSNames = pathNames
		numFZ = len(featZSNames)
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
		gKnown = mp.readFileAsList(si + 'known.txt')
		giKnown = mp.convertToIndices(gKnown, geneDict)
		sumFTV = np.sum(featTWVals[giKnown,:], axis=0)
		keepIdx = np.nonzero(sumFTV)
#		print("\n {} \n".format(keepIdx[0]))
		numFT = len(keepIdx[0])

		print("    ... including term membership features")
		features = np.hstack( (features, featTWVals[:,keepIdx[0]]) )
		featNames.extend(np.ravel(featTWNames[keepIdx]))
	# verify some features have been loaded
	if features.shape[1] == 0 :
		print("ERROR: No features were specified for classification.")
		sys.exit 
	#end if
	numFAll = numFP + numFZ + numFN + numFT


	# Normalize the feature values
	features = mp.normalizeFeatureColumns(features)



	# 6) Prepare the test/train vectors & labels
	print("  Creating train & test data ...")

	# Cluster to create train/test sets

	clusLabel, featLabel = mp.clusterTrainSets(si, geneDict, features, maxClusters)
	clusVals = np.unique(clusLabel)	
	nClus = len(clusVals) - 1	# less 1 b/c first label is the Known set
	if newVerbose :
		print("  The Unknown set was grouped into {} clusters.".format(nClus))
	#end if

	# For each cluster, use label 0 + label N to train
	#	then predict on the rest
	geneRanks = np.zeros( (len(geneDict), nClus), dtype=np.int32 )
	geneScores = np.zeros( (len(geneDict), nClus), dtype=np.float32 )

	# Store the number of clusters which use certain features
	featT1Dict = dict()
	featT1Set = set()
	featT5Dict = dict()
	featT5Set = set()
	featTADict = dict()
	featTASet = set()

	col = -1
	for cv in clusVals :
		# Skip cv == 0 (will use 0 as pos train set)
		if cv == 0 :
			continue
		col += 1

		# Extract the train set & labels
		trainIdx = np.hstack(( np.where(clusLabel == 0)[0], np.where(clusLabel == cv)[0] ))
		trainIdx.sort()
		trainSet = features[trainIdx,:]
		trainLabel = featLabel[trainIdx]
		trainLabel = np.ravel(trainLabel)

		# Extract the test set & labels
		testIdx = [idx for idx in range(len(clusLabel)) if idx not in trainIdx]
		testIdx = np.asarray(testIdx)
		# error hack when testing code with small samples :
		if len(testIdx) == 0 :
			testIdx = [0]
		#end if
		testSet = features[testIdx,:]
		testLabel = featLabel[testIdx]
		testLabel = np.ravel(testLabel)
	#end loop



# TODO: If the classifier has 0 coefficients, re-run? discount?

		# 7) Train the classifier
		numCoefs = 0
		loopCount = 0
		while numCoefs == 0:
			if useCfier == 1 :	# 1 = Lasso
				cfier = lm.LassoCV(alphas=useGivenRange, positive=usePos,
					max_iter=lMaxIter, normalize=lNorm, fit_intercept=lFitIcpt)
			elif useCfier == 2 :	# 2 = ElasticNet
				cfier = lm.ElasticNetCV(l1_ratio=enRatios, positive=usePos, fit_intercept=enFitIncept,
					n_alphas=enNAlphas, normalize=enNorm, copy_X=enCopy, max_iter=enMaxIter)
			elif useCfier == 3 :	#3 = SVM
				cfier = svm.SVC(probability=True, max_iter=svmMaxIter, decision_function_shape='ovr')
			elif useCfier == 4 :	#4 = LinearSVM
				cfier = svm.LinearSVC(penalty='l1', max_iter=svmMaxIter, dual=False)
			elif useCfier == 5 :	# 5 = Lasso + SVM
				# part 1: get coefficients
				cfier = lm.LassoCV(alphas=useGivenRange, positive=usePos,
					max_iter=lMaxIter, normalize=lNorm, fit_intercept=lFitIcpt)
			else :
				print("ERROR: specified classifier unrecognized: useCfier = {}".format(useCfier))
			#end if

			cfier.fit(trainSet, trainLabel)
			if useCfier == 3 :
#				print(cfier.support_vectors_)
#				print(cfier.n_support_)
				numCoefs = numFAll
			else :
				numCoefs = len(np.nonzero(cfier.coef_)[0])
			#end if

			loopCount += 1
			if loopCount >= 3 :
				break
		#end loop

		if newVerbose :
			# view quick statistics from this training session
			if useCfier < 3 :
				print("Clus {}; cfier {}; pos={}; score: {:.3f}; coeffs: {}".format(int(cv), useCfier,
					usePos, cfier.score(trainSet, trainLabel), len(np.nonzero(cfier.coef_)[0])))
				print("  iterations: {}; chosen alpha: {:.6f}".format(cfier.n_iter_, cfier.alpha_))
				if useCfier == 2 :	# 2 = ElasticNet
					print("    l1 ratio: {}".format( cfier.l1_ratio_ ))
			elif useCfier == 3 :
#TODO: this
				print("Clus {}; cfier {}; score: {:.3f}; coeffs: {}".format(int(cv), useCfier,
					cfier.score(trainSet, trainLabel), numFAll))
#				print("  iterations: {}".format(cfier.n_iter_))
			elif useCfier == 4 :
					print("Clus {}; cfier {}; score: {:.3f}; coeffs: {}".format(int(cv), useCfier,
					cfier.score(trainSet, trainLabel), len(np.nonzero(cfier.coef_)[0])))
			elif useCfier == 5 :
				useFeatIdx = np.nonzero(cfier.coef_)[0]
				print("Clus {}; cfier {}; pos={}; score: {:.3f}; coeffs: {}".format(int(cv), useCfier,
					usePos, cfier.score(trainSet, trainLabel), len(useFeatIdx)))
				print("  iterations: {}; chosen alpha: {:.6f}".format(cfier.n_iter_, cfier.alpha_))

			#end if
			if loopCount > 1 :
				print("    attempts: {}".format(loopCount))
		#end if



		# 8) Predict on the test set
		#	score the genes, then sort, place rank into array geneRanks
		ranker = np.recarray(testSet.shape[0],
			dtype=[('inverse', 'f4'), ('score', 'f4'), ('nameIdx', 'a20')])

		if useCfier < 3 :
			cfPredLabel = cfier.predict(testSet)
		elif useCfier == 3 :
#			cfPredLabel = cfier.predict(testSet)
			cfPredLabel = cfier.decision_function(testSet)
		elif useCfier == 4 :
			cfPredLabel = cfier.decision_function(testSet)
		elif useCfier == 5 :
			if numCoefs > 0 :
				cfier = svm.LinearSVC(penalty='l2', max_iter=svmMaxIter, dual=False)
				cfier.fit(trainSet[:,useFeatIdx], trainLabel)
				cfPredLabel = cfier.decision_function(testSet[:,useFeatIdx])
		#end if
		cfPredLabel = np.ravel(cfPredLabel)

		ranker['score'] = cfPredLabel
		ranker['inverse'] = 0 - cfPredLabel
		ranker['nameIdx'] = testIdx

		ranker.sort(order=['inverse', 'nameIdx'])

		# If no coeffs were chosen, rank will be random; leave col as all zeros
		if numCoefs > 0 :
			# Place genes' rank & score into appropriate matrices
			rank = 0
			for entry in ranker :
				rank += 1
				row = int(entry['nameIdx'])
				geneRanks[row, col] = rank
				geneScores[row,col] = entry['score']
			#end loop
		#end if



		# 9) Collect top features (& other stats?)

		# For each cluster, for each feature,
		#	increment count if feat in top X for that cluster
		#	will later output how often feat used across all clusters

		# Extract indices corresponding to top 5 weighted features
		if useCfier != 3 :

			if useCfier > 3 :
				featWeights = cfier.coef_[0]
			else :
				featWeights = cfier.coef_
#			print(cfier.coef_)
#			print(np.nonzero(featWeights))
			numFeats = len(np.nonzero(featWeights)[0])
			if numFeats > 0 :
				topFeats = np.ones( (numFeats), dtype=np.int32 ) * (-1)
				for num in range(numFeats) :
					featIdx = np.argmax(featWeights)
					topFeats[num] = featIdx
					featWeights[featIdx] = 0
				#end loop

				# Increment count for the Top 1 path
				if topFeats[0] in featT1Set :
					featT1Dict[topFeats[0]] += 1
				else :
					featT1Dict[topFeats[0]] = 1
					featT1Set.add(topFeats[0])
				#end if

				# Increment count for the Top 5 paths
				for num in range(5) :
					if numFeats <= num :
						break
					#end if
					if topFeats[num] in featT5Set :
						featT5Dict[topFeats[num]] += 1
					else :
						featT5Dict[topFeats[num]] = 1
						featT5Set.add(topFeats[num])
				#end loop

				# Increment count for all non-zero paths
				for num in range(numFeats) :
					if numFeats <= num :
						break
					#end if
					if topFeats[num] in featTASet :
						featTADict[topFeats[num]] += 1
					else :
						featTADict[topFeats[num]] = 1
						featTASet.add(topFeats[num])
				#end loop
			#end if
		#end if

	#end loop (per-cluster loop)



	# 10) Output the ranked genes to file

	# Get the AVERAGE gene rank
	# First, vote across the clusters and sort by rank
	# NOTE: any genes in the Known Set should have rank == 0
	cfGenes = np.recarray( len(geneDict), dtype=[('geneIdx', 'i4'), ('rank', 'f4')])
	cfGenes['geneIdx'] = np.arange(len(geneDict))
	cfGenes['rank'] = np.mean(geneRanks, axis=1)
	cfGenes.sort(order=['rank','geneIdx'])

	# write the file
	fname = 'ranked_genes-' + useLabel + '_Avg.txt'
	print("  Saving ranked genes to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		firstRow = True
		for row in range(len(cfGenes)) :
			# Skip any genes whose rank == 0 (the Known Set)
			if cfGenes['rank'][row] == 0 :
				continue
			#end if
			if not firstRow :
				fout.write('\n')
			fout.write('{:3.3f}{}{}'.format(cfGenes['rank'][row],
				textDelim, geneNames[cfGenes['geneIdx'][row]]))
			firstRow = False
	#end with

	# Get the MEDIAN gene rank
	# First, vote across the clusters and sort by rank
	# NOTE: any genes in the Known Set should have rank == 0
	cfGenes = np.recarray( len(geneDict), dtype=[('geneIdx', 'i4'), ('rank', 'f4')])
	cfGenes['geneIdx'] = np.arange(len(geneDict))
	cfGenes['rank'] = np.median(geneRanks, axis=1)
	cfGenes.sort(order=['rank','geneIdx'])

	# write the file
	fname = 'ranked_genes-' + useLabel + '_Med.txt'
	print("  Saving ranked genes to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		firstRow = True
		for row in range(len(cfGenes)) :
			# Skip any genes whose rank == 0 (the Known Set)
			if cfGenes['rank'][row] == 0 :
				continue
			#end if
			if not firstRow :
				fout.write('\n')
			fout.write('{:3.3f}{}{}'.format(cfGenes['rank'][row],
				textDelim, geneNames[cfGenes['geneIdx'][row]]))
			firstRow = False
	#end with


	# 11) Output the selected feature info to file
	if useCfier != 3 :

		# Sort the Top 1 paths
		featT1Sort = np.recarray( len(featT1Dict), dtype=[('pathIdx', 'i4'), ('count', 'i4')])
		pathIdxList = list(featT1Dict.keys())
		row = -1
		for item in pathIdxList :
			row += 1
			featT1Sort['pathIdx'][row] = item
			featT1Sort['count'][row] = featT1Dict[item]
		#end if
		featT1Sort[::-1].sort(order=['count', 'pathIdx'])

		# Save the Top 1 paths to file
		fname = 'ranked_features_Top1-' + useLabel + '.txt'
		with open(si + fname, 'w') as fout :
			fout.write('Clusters:{}{}'.format(textDelim, nClus))
			for row in range(len(featT1Sort)) :
				fout.write('\n{}{}{}'.format(featT1Sort['count'][row],
					textDelim, featNames[featT1Sort['pathIdx'][row]]))
		#end with

		# Sort the Top 5 paths
		featT5Sort = np.recarray( len(featT5Dict), dtype=[('pathIdx', 'i4'), ('count', 'i4')])
		pathIdxList = list(featT5Dict.keys())
		row = -1
		for item in pathIdxList :
			row += 1
			featT5Sort['pathIdx'][row] = item
			featT5Sort['count'][row] = featT5Dict[item]
		#end if
		featT5Sort[::-1].sort(order=['count', 'pathIdx'])

		# Save the Top 5 paths to file
		fname = 'ranked_features_Top5-' + useLabel + '.txt'
		with open(si + fname, 'w') as fout :
			fout.write('Clusters:{}{}'.format(textDelim, nClus))
			for row in range(len(featT5Sort)) :
				fout.write('\n{}{}{}'.format(featT5Sort['count'][row],
					textDelim, featNames[featT5Sort['pathIdx'][row]]))
		#end with

		# Sort the Top All Non-Zero paths
		featTASort = np.recarray( len(featTADict), dtype=[('pathIdx', 'i4'), ('count', 'i4')])
		pathIdxList = list(featTADict.keys())
		row = -1
		for item in pathIdxList :
			row += 1
			featTASort['pathIdx'][row] = item
			featTASort['count'][row] = featTADict[item]
		#end if
		featTASort[::-1].sort(order=['count', 'pathIdx'])

		# Save the Top All Non-Zero paths to file
		fname = 'ranked_features_TopNZ-' + useLabel + '.txt'
		with open(si + fname, 'w') as fout :
			fout.write('Clusters:{}{}'.format(textDelim, nClus))
			for row in range(len(featTASort)) :
				fout.write('\n{}{}{}'.format(featTASort['count'][row],
					textDelim, featNames[featTASort['pathIdx'][row]]))
		#end with
	#end if


	# 12) Output the parameters to file
	fname = 'parameters-' + useLabel + '.txt'
	with open(si+fname, 'w') as fout :
		fout.write('\n')
		fout.write('Sampling Method for Neg examples\n')
		fout.write('  as One-Class w/ clustering on the Unknown set\n')
		fout.write('\n')
		fout.write('Features Used\n')
		fout.write('PathSim sum:{}{}\n'.format(textDelim, useFeatPathSim))
		fout.write('path Z-Score:{}{}\n'.format(textDelim, useFeatPathZScore))
		fout.write('Neighborhood:{}{}\n'.format(textDelim, useFeatNeighbor))
		fout.write('Term Weights:{}{}\n'.format(textDelim, useFeatTermWeights))
		fout.write('\n')

#TODO: collect some stats (ie: common alphas, l1 ratios, etc)
		fout.write('Classifier Parameters\n')
		if useCfier == 1 :
			fout.write('method:{}Lasso\n'.format(textDelim))
			fout.write('positive:{}{}\n'.format(textDelim, usePos))
			# fout.write('alpha range:{}{}\n'.format(textDelim, useGivenRange))
			# fout.write('alpha chosen:{}{}\n'.format(textDelim, cfier.alpha_))
			fout.write('max_iter:{}{}\n'.format(textDelim, lMaxIter))
			fout.write('normalize:{}{}\n'.format(textDelim, lNorm))
			fout.write('fit_intercept:{}{}\n'.format(textDelim, lFitIcpt))
		elif useCfier == 2 :
			fout.write('method:{}ElasticNet\n'.format(textDelim))
			fout.write('ratio range:{}{}\n'.format(textDelim, enRatios))
			# fout.write('ratio chosen:{}{}\n'.format(textDelim, cfENet.l1_ratio_))
			fout.write('num alphas:{}{}'.format(textDelim, enNAlphas))
			fout.write('max_iter:{}{}\n'.format(textDelim, enMaxIter))
			fout.write('normalize:{}{}\n'.format(textDelim, enNorm))
			fout.write('fit_intercept:{}{}\n'.format(textDelim, enFitIncept))
		elif useCfier == 3 :
			fout.write('method:{}SVM\n'.format(textDelim))
			fout.write('kernel:{}{}\n'.format(textDelim, svmKernel))
			fout.write('max_iter:{}{}\n'.format(textDelim, svmMaxIter))
		fout.write('\n')
#TODO: output for cfiers 4 & 5
		fout.write('Prediction Results\n')
		if useCfier == 3 :
			fout.write('Training score:{}{:3.3f}\n'.format(textDelim, cfier.score(trainSet, trainLabel)))
			fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, cfier.score(testSet, testLabel)))
		elif useCfier == 5 :
			fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(useFeatIdx)))
			fout.write('Training score:{}{:3.3f}\n'.format(textDelim, cfier.score(trainSet[:,useFeatIdx], trainLabel)))
			fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, cfier.score(testSet[:,useFeatIdx], testLabel)))
		else :
			fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(np.nonzero(cfier.coef_)[0])))
			fout.write('Training score:{}{:3.3f}\n'.format(textDelim, cfier.score(trainSet, trainLabel)))
			fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, cfier.score(testSet, testLabel)))
		fout.write('\n')
	#end with


	thisRound += 1
	print("--{} of {}".format(thisRound, len(dSubDirs)))
	print("--elapsed time: {:.3} (s)".format(time.time()-tstart))
#end loop



print("\nDone.\n")