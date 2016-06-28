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
#	along side the PathSim ones
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
sDir = '../Dropbox/mp/output/pred04-set01'

# File name containing feature vectors
fSimilarity = 'features_PathSim.gz' 

#retCutoffs = [50, 100, 200, 500, 1000, 2000]


# Adjustible classifier parameters
useCfier = 1
	# 1 = Lasso, 2 = ElasticNet, ...
usePos = True
	# True/False: limit to only Positive coefficient values
useFeatPaths = True
	# True/False: use the pathsim sum features
useFeatNeighbor = True
	# True/False: use the neighborhood features
useGivenRange = np.linspace(0.00005, 0.002, num=13)
	# array of vals; 'None' means to auto-search for alphas
# useLabel = 'ClusVote_cLasP_fPN'
# 	# string: label for the output files
# 	# ie: ClusVote_c<Las/Enet/Log><P for Pos>_f<P for pathsim><N for neighborhood>
# #TODO: automate useLabel based on params


# LASSO params
lMaxIter = 2000
lNorm = True
lFitIcpt = True

# Elastic Net params
enRatios = [0.2, 0.4, 0.6, 0.75, 0.85, 0.9, 0.95]
enNAlphas = 17
enMaxIter = 5000
enFitIncept = True
enNorm = True
enCopy = True


# verbose feedback ?
newVerbose = True

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
useLabel = 'ClusVote_c'
if useCfier == 1 :
	useLabel = useLabel + 'Las'
elif useCfier == 2 :
	useLabel = useLabel + 'Enet'
else :
	print("ERROR: useCfier value is unrecognized: {}".format(useCfier))
if usePos :
	useLabel = useLabel + 'Pos'
useLabel = useLabel + '_f'
if useFeatPaths :
	useLabel = useLabel + 'P'
if useFeatNeighbor :
	useLabel = useLabel + 'N'
#end if



# 1) Load the gene-index dictionary & path names
geneDict, pathDict = mp.getGeneAndPathDict(sDir)
geneNames = list(geneDict.keys())
geneNames.sort()
pathNames = mp.removeInvertedPaths(pathDict)
del pathDict



# 2) Load the neighborhood features
if useFeatNeighbor :
	featNbVals, featNbNames = mp.getFeaturesNeighborhood(sDir, 'LogScale')
	featNbNames = np.ravel(featNbNames)
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
	if useFeatPaths :
		featPSVals = mp.readFileAsMatrix(si, fSimilarity)
		# NOTE: previous version of mpPredict04a had extra zeros at end of vectors
		#	discarding those columns
		featPSVals = featPSVals[:,0:len(pathNames)]
		featPSNames = pathNames
#		featPSNames = np.reshape(featPSNames, (1, len(featPSNames)))
	#end if



	# 5) Combine the features as specified by parameters (useFeat...)
	features = np.zeros( (len(geneDict), 0), dtype=np.float32)
	featNames = list()

#	print(featNames)

#	print("{},{}".format(features.shape, featPSVals.shape))
	if useFeatPaths :
		print("    ... including PathSim sum features")
		features = np.hstack( (features, featPSVals) )
		featNames.extend(featPSNames)
#		featNames.extend(np.ravel(featPSNames))
#		featNames = np.hstack( (featNames, featPSNames))
#		print(len(featPSNames))
#		print(len(featNames))
	if useFeatNeighbor :
		print("    ... including neighborhood features")
		features = np.hstack( (features, featNbVals) )
		featNames.extend(np.ravel(featNbNames))
#		featNames = np.hstack( (featNames, featNbNames) )
#		print(len(featNbNames))
#		print(len(featNames))
	# verify some features have been loaded
	if features.shape[1] == 0 :
		print("ERROR: No features were specified for classification.")
		sys.exit 
	#end if



	# 6) Prepare the test/train vectors & labels
	print("  Creating train & test data ...")

# 	if useFeatNeighbor :
# 		print("    ... including neighborhood features")
# 		features = np.hstack( (featNbVals, featPSVals) )
# 		featNames = np.hstack( (featNbNames, featPSNames) )
# 	else :
# #TODO: make pathsim features optional
# 		features = featPSVals
# 		featNames = featPSNames
# 	#end if

	# Cluster to create train/test sets
#	trainSet, trainLabel, testSet, testLabel, giTest = mp.createTrainTestSets(si, geneDict, features, True)

#	print(np.amax(features[:,0]))
	clusLabel, featLabel = mp.clusterTrainSets(si, geneDict, features)
#	print(np.amax(features[:,0]))
#	print(np.unique(featLabel))

	clusVals = np.unique(clusLabel)	
#	nClus = np.amax(clusLabel)
	nClus = len(clusVals) - 1	# less 1 b/c first label is the Known set
	if newVerbose :
#		print("verify: amax={}, (unique-1)={}".format(
#			np.amax(clusLabel), len(clusVals) - 1))
		print("  The Unknown set was grouped into {} clusters.".format(nClus))
#	print(np.unique(clusLabel))


	# For each cluster, use label 0 + label N to train
	#	then predict on the rest
	geneRanks = np.zeros( (len(geneDict), nClus), dtype=np.int32 )
	geneScores = np.zeros( (len(geneDict), nClus), dtype=np.float32 )
	featT1Dict = dict()
	featT1Set = set()
	featT5Dict = dict()
	featT5Set = set()

	col = -1
	for cv in clusVals :
		# Skip cv == 0 (will use 0 as pos train set)
		if cv == 0 :
			continue

		col += 1

		# Extract the train set & labels
		trainIdx = np.hstack(( np.where(clusLabel == 0)[0], np.where(clusLabel == cv)[0] ))
		trainIdx.sort()
#		print(trainIdx.shape)
		trainSet = features[trainIdx,:]
		trainLabel = featLabel[trainIdx]
		trainLabel = np.ravel(trainLabel)

		# Extract the test set & labels
		testIdx = [idx for idx in range(len(clusLabel)) if idx not in trainIdx]
		testIdx = np.asarray(testIdx)
#		testIdx = np.reshape(testIdx, (len(testIdx),1))
#		print(testIdx.shape)
		testSet = features[testIdx,:]
		testLabel = featLabel[testIdx]
		testLabel = np.ravel(testLabel)
#		testNames = geneNames[testIdx]



		# 7) Train the classifier
		if useCfier == 1 :	# 1 = Lasso
			cfier = lm.LassoCV(alphas=useGivenRange, positive=usePos,
				max_iter=lMaxIter, normalize=lNorm, fit_intercept=lFitIcpt)
			# cfier = lm.LassoCV(positive=usePos,
			# 	max_iter=lMaxIter, normalize=lNorm, fit_intercept=lFitIcpt)
		elif useCfier == 2 :	# 2 = ElasticNet
#TODO: this
			cfier = lm.ElasticNetCV()
		else :
			print("ERROR: specified classifier unrecognized: useCfier = {}".format(useCfier))
		#end if

		cfier.fit(trainSet, trainLabel)

		if newVerbose :
			# view quick statistics from this training session
			print("Clus {}; cfier {}; pos={}; score: {:.3f}; coeffs: {}".format(cv, useCfier,
				usePos, cfier.score(trainSet, trainLabel), len(np.nonzero(cfier.coef_)[0])))
			print("  iterations: {}; chosen alpha: {:.6f}".format(cfier.n_iter_, cfier.alpha_))
			if useCfier == 2 :	# 2 = ElasticNet
				print("    l1 ratio: {}".format( cfier.l1_ratio_ ))
		#end if



		# 8) Predict on the test set
		#	score the genes, then sort, place rank into array geneRanks
		ranker = np.recarray(testSet.shape[0],
#			dtype=[('inverse', 'f4'), ('score', 'f4'), ('names', 'a20')])
			dtype=[('inverse', 'f4'), ('score', 'f4'), ('nameIdx', 'a20')])


		cfPredLabel = cfier.predict(testSet)
		cfPredLabel = np.ravel(cfPredLabel)

		ranker['score'] = cfPredLabel
		ranker['inverse'] = 0 - cfPredLabel
#		ranker['names'] = testNames
		ranker['nameIdx'] = testIdx
#		for r in len(ranker) :
#			row = [ 0 - cfPredLabel[r], cfPredLabel[r], r]
#			ranker[r] = row
#		#end loop

#		ranker.sort(order=['inverse', 'names'])
		ranker.sort(order=['inverse', 'nameIdx'])

		rank = 0
		for entry in ranker :
			rank += 1
#			row = geneDict[entry['names']]
			row = int(entry['nameIdx'])
			geneRanks[row, col] = rank
			geneScores[row,col] = entry['score']
		#end loop

#		print(geneRanks[row,col])


		# 9) 
#TODO: collect feature names & meaningful stats

		# Extract indices corresponding to top 5 weighted features
		featWeights = cfier.coef_
		topPaths = np.ones( (5), dtype=np.int32 ) * (-1)
		for num in range(5) :
			if (len(np.nonzero(featWeights)[0]) <= num) :
				break
			featIdx = np.argmax(featWeights)
			topPaths[num] = featIdx
			featWeights[featIdx] = 0
		#end if

		# Increment count for the Top 1 path
		if topPaths[0] in featT1Set :
			featT1Dict[topPaths[0]] += 1
		else :
			featT1Dict[topPaths[0]] = 1
			featT1Set.add(topPaths[0])
		#end if

		# Increment count for the Top 5 paths
		for num in range(5) :
			if topPaths[num] in featT5Set :
				featT5Dict[topPaths[num]] += 1
			else :
				featT5Dict[topPaths[num]] = 1
				featT5Set.add(topPaths[num])
		#end loop


#		break
	#end loop


#	break



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
#TODO: this!

	# Sort the Top 1 paths
	featT1Sort = np.recarray( len(featT1Dict), dtype=[('pathIdx', 'i4'), ('count', 'i4')])
	pathIdxList = list(featT1Dict.keys())
	row = -1
	for item in pathIdxList :
		row += 1
		featT1Sort['pathIdx'] = item
		featT1Sort['count'] = featT1Dict[item]
	#end if
	featT1Sort[::-1].sort(order=['count', 'pathIdx'])

#	print(featT1Sort)
#	print(len(featNames))
#	print(featNames)

	# Save the Top 1 paths to file
	fname = 'ranked_features_Top1-' + useLabel + '.txt'
	with open(si + fname, 'w') as fout :
#		fout.write('intercept:{}{}'.format(textDelim, cfier.intercept_))
		for row in range(len(featT1Sort)) :

#			print(featT1Sort['pathIdx'][row])
#			print(featT1Sort[row]['pathIdx'])
#			print('{},{}'.format(row, featT1Sort['pathIdx'][row]))
#			print('{},{}'.format(row, featT1Sort['count'][row]))

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
#		fout.write('intercept:{}{}'.format(textDelim, cfier.intercept_))
		for row in range(len(featT5Sort)) :
			fout.write('\n{}{}{}'.format(featT5Sort['count'][row],
				textDelim, featNames[featT5Sort['pathIdx'][row]]))
	#end with



	# 12) Output the parameters to file
	fname = 'parameters-' + useLabel + '.txt'
	with open(si+fname, 'w') as fout :
		fout.write('\n')
		fout.write('Sampling Method for Neg examples\n')
		fout.write('  as One-Class w/ clustering on the Unknown set\n')
		fout.write('\n')
		fout.write('Features Used\n')
		fout.write('PathSim sum:{}{}\n'.format(textDelim, useFeatPaths))
		fout.write('Neighborhood:{}{}\n'.format(textDelim, useFeatNeighbor))
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
		fout.write('\n')

#		fout.write('Similarity Metric:{}PathSim sum over set\n'.format(textDelim))
		fout.write('Prediction Results\n')
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