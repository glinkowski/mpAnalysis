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
sDir = '../Dropbox/mp/output/pred04-test01'

# File name containing feature vectors
fSimilarity = 'features_PathSim.gz' 

retCutoffs = [50, 100, 200, 500, 1000, 2000]


useCfier = 1
	# 1 = Lasso, 2 = ElasticNet, ...
usePos = True
	# True/False
useFeatPaths = True
	# True/False: use the pathsim sum features
useFeatNeighbor = True
	# True/False: use the neighborhood features
useGivenRange = np.linspace(0.00005, 0.002, num=23)
	# array of vals; 'None' means to auto-search for alphas
useLabel = 'Lasso_Cluster_cP_fPN'
	# string: label for the output files


# LASSO params
#lAlpha02 = [0.0008, 0.0007, 0.0006, 0.0005, 0.0004, 0.0003]
#lAlphaPos = [0.001, 0.0007, 0.0006, 0.0005, 0.0004, 0.0003, 0.0001, 0.00008]
lMaxIter = 3000
lNorm = True
lFitIcpt = True


# Elastic Net params
#enRatios = [0.4, 0.7, 0.85, 0.9, 0.95]
enRatios = [0.2, 0.4, 0.6, 0.75, 0.85, 0.9, 0.95]
enNAlphas = 17
enMaxIter = 5000
enFitIncept = True
enNorm = True
enCopy = True


# # Lables for pos & neg training sets
# pLabel = 1
# nLabel = 0

# verbose feedback ?
newVerbose = True

textDelim = '\t'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()
print("\nPerforming regression(s) on {}".format(sDir))

mp.setParamVerbose(newVerbose)



# 1) Load the gene-index dictionary & path names
geneDict, pathDict = mp.getGeneAndPathDict(sDir)
geneNames = list(geneDict.keys())
geneNames.sort()
pathNames = mp.removeInvertedPaths(pathDict)
del pathDict


# 2) Load the neighborhood features
featNbVals, featNbNames = mp.getFeaturesNeighborhood(sDir, 'LogScale')


# 3) Loop over all the subdirectories
dSubDirs = mp.getSubDirectoryList(sDir)

thisRound = 0
#for si in dSubDirs[0:1] :
for si in dSubDirs :

	# Display directory to examine
	sv = si.split('/')
	print("\n{}/{}/".format(sv[-3],sv[-2]))


	# 4) Load the PathSim features
	featPSVals = mp.readFileAsMatrix(si, fSimilarity)
	# NOTE: previous version of mpPredict04a had extra zeros at end of vectors
	#	discarding those columns
	featPSVals = featPSVals[:,0:len(pathNames)]
	featPSNames = pathNames
	featPSNames = np.reshape(featPSNames, (1, len(featPSNames)))



	# 5) Prepare the test/train vectors & labels
	print("  Creating train & test data ...")

	if useFeatNeighbor :
		print("    ... including neighborhood features")
		features = np.hstack( (featNbVals, featPSVals) )
		featNames = np.hstack( (featNbNames, featPSNames) )
	else :
#TODO: make pathsim features optional
		features = featPSVals
		featNames = featPSNames
	#end if

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
		print("verify: amax={}, (unique-1)={}".format(
			np.amax(clusLabel), len(clusVals) - 1))
	print("  The Unknown set was grouped into {} clusters.".format(nClus))
#	print(np.unique(clusLabel))


	# For each cluster, use label 0 + label N to train
	#	then predict on the rest
	geneRanks = np.zeros( (len(geneDict), nClus), dtype=np.int32 )

	col = -1
	for cv in clusVals :
		col += 1

		# Skip cv == 0 (will use 0 as pos train set)
		if cv == 0 :
			continue

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
#		testLabel = featLabel[testIdx]
#		testLabel = np.ravel(testLabel)
#		testNames = geneNames[testIdx]

		# 6) Train the classifier
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

#TODO: collect feature names & stats

		# 7) Predict on the test set
		#	score the genes, then sort, place rank into array geneRanks
		geneScores = np.recarray(testSet.shape[0],
#			dtype=[('inverse', 'f4'), ('score', 'f4'), ('names', 'a20')])
			dtype=[('inverse', 'f4'), ('score', 'f4'), ('nameIdx', 'a20')])


		cfPredLabel = cfier.predict(testSet)
		cfPredLabel = np.ravel(cfPredLabel)

		geneScores['score'] = cfPredLabel
		geneScores['inverse'] = 0 - cfPredLabel
#		geneScores['names'] = testNames
		geneScores['nameIdx'] = testIdx
#		for r in len(geneScores) :
#			row = [ 0 - cfPredLabel[r], cfPredLabel[r], r]
#			geneScores[r] = row
#		#end loop

#		geneScores.sort(order=['inverse', 'names'])
		geneScores.sort(order=['inverse', 'nameIdx'])

		rank = 0
		for entry in geneScores :
			rank += 1
#			row = geneDict[entry['names']]
			row = int(entry['nameIdx'])
			geneRanks[row, col] = rank
		#end loop

		print(geneRanks[row,col])


		break
	#end loop


	break

#TODO: 
#	cluster the Unknown, get those labels
#	for each cluster, create train/test sets
#	train the classifier
#	return test prediction results
#	Then: vote, collect & output results






	# 5) Output results to file, Lasso 2-class

	# Save the selected paths & scores/weights
	# 	feature coefficients are the metapath weights
	cfCoefs = np.nonzero(cfLassoStd.coef_)[0]
	cfPaths = np.recarray( len(cfCoefs), dtype=[('path', 'i4'), ('weight', 'f4')] )
	for row in range(len(cfCoefs)) :
		cfPaths[row] = (row, cfCoefs[row])
	cfPaths[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

	# write the file
	fname = 'ranked_paths-Lasso_2ClassStd.txt'
	print("Saving data for the Lasso standard approach ...")
	print("  Saving top paths to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		fout.write('intercept:{}{}'.format(textDelim, cfLassoStd.intercept_))
		for row in range(len(cfPaths)) :
			fout.write('\n{}{}{}'.format(cfPaths['weight'][row],
				textDelim, pathNames[cfPaths['path'][row]]))
	#end with

	#Sort the genes by (inverse) rank
	cfGenes = np.recarray( len(cfPredLabel), dtype=[('gene', 'i4'), ('rank', 'f4')] )
	cfGenes['gene'] = giTest
	cfGenes['rank'] = cfPredLabel
	cfGenes[::-1].sort(order=['rank','gene'])

	# write the file
	fname = 'ranked_genes-Lasso_2ClassStd.txt'
	print("  Saving ranked genes to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		firstRow = True
		for row in range(len(cfGenes)) :
			if not firstRow :
				fout.write('\n')
			fout.write('{:3.3f}{}{}'.format(cfGenes['rank'][row],
				textDelim, geneNames[cfGenes['gene'][row]]))
			firstRow = False
	#end with

	# Save the parameters & results
	fname = 'parameters-Lasso_2ClassStd.txt'
	with open(si+fname, 'w') as fout :
		fout.write('\n')
		fout.write('Sampling Method for Neg examples\n')
		fout.write('  as Two-Class\n')
		fout.write('\n')

		fout.write('Lasso Parameters\n')
		fout.write('method:{}Lasso (standard)\n'.format(textDelim))
		fout.write('a_range:{}{}\n'.format(textDelim, lAlpha02))
		fout.write('alpha:{}{}\n'.format(textDelim, cfLassoStd.alpha_))
		fout.write('max_iter:{}{}\n'.format(textDelim, lMaxIter))
		fout.write('normalize:{}{}\n'.format(textDelim, lNorm))
		fout.write('positive:{}{}\n'.format(textDelim, 'False'))
		fout.write('fit_intercept:{}{}\n'.format(textDelim, lFitIcpt))
		fout.write('\n')

		fout.write('Similarity Metric:{}PathSim sum over set\n'.format(textDelim))
		fout.write('Prediction Results\n')
		fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(np.nonzero(cfLassoStd.coef_)[0])))
		fout.write('Training score:{}{:3.3f}\n'.format(textDelim, cfLassoStd.score(trainSet, trainLabel)))
		fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, cfLassoStd.score(testSet, testLabel)))
		fout.write('\n')
	#end with



	thisRound += 1
	print("    --{} of {}".format(thisRound, len(dSubDirs)))
	print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))
#end loop



print("\nDone.\n")