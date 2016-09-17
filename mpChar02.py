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

import mpCharLib as cl
import time
import numpy as np
import gzip
from sklearn import linear_model as lm
import random



######## ######## ######## ######## 
# PARAMETERS

# folder containing the pre-processed samples
sDir = '../Dropbox/mp/output/pred04-dbgap300'

# Control the iterations & error
numIterations = 1
	# how many iterations to perform
numVotes = 11
	# how many random samples for comparison
# numWrong = 3
# 	# how many False Negatives before next iteration
numExitKnown = 29
	# minimum number of Known genes before iterations stop
numExitUnknown = 399
	# min number of Unknonw genes before iterations halt

retryOnZeroCoeffs = True
	# whether to allow coeffs = 0 in results
retrySubPortion = 0.75
	# how many of Known to keep in new sub-sample
retryMinValid = 9
	# minimum Known genes to use for PosTrain

# verbose feedback ?
verbose = True



#TODO: remove / hardcode the following parameters


# Adjustible classifier parameters
useCfier = 2
	# 1 = Lasso, 2 = ElasticNet, ...
limitMPLen = [1, 2, 3]
	# select only meta-paths of specific length(s)
usePos = True
	# True/False: limit to only Positive coefficient values
useFeatPathSim = False
	# True/False: use the pathsim sum features
fSimilarity = 'features_PathSim.gz'
	# File name containing pathsim feature vectors
useFeatPathZScore = True
	# True/False: use the pathsim sum features
fZScoreSim = 'features_ZScoreSim.gz'
	# File name containing path z-score vectors
useFeatTermWeights = False
	# True/False: use the indirect term features
useFeatNeighbor = False
	# True/False: use the neighborhood features
useGivenRange = np.linspace(0.00001, 0.05, num=27)
	# array of vals; 'None' means to auto-search for alphas
# maxClusters = 11
# 	# maximum number of clusters to use
# 	#	-1 = no maximum value is set

# Label for pos & neg labels
pLabel = 1
nLabel = 0
negMultiplier = 1

# LASSO params
lMaxIter = 800
lNorm = True
lFitIcpt = True

# Elastic Net params
enRatios = [0.2, 0.4, 0.6, 0.75, 0.85, 0.9, 0.95]
enNAlphas = 11
enMaxIter = 400
enFitIncept = True
enNorm = True
enCopy = True



# verbose feedback ?
verboseOutput = True

######## ######## ######## ######## 




######## ######## ######## ######## 
# PRIMARY FUNCTION
#TODO: add description
def predictIterative(printFlag) :

	# Parameters
	numVotes = 11
	usePos = True
	negMultiplier = 1


	if printFlag :
		print("\nPerforming regression(s) on {}".format(sDir))



	# 0) Create the useLabel variable
	# string: label for the output files
	# ie: ClusVote_c<Las/Enet/Log/SVM><P for Pos>_f<P for pathsim><Z for z-score>
	#	<T for term weights><N for neighborhood>_m<neg sample size multiplier>
	useLabel = 'Vote{}_cLas'.format(numVotes)
	if usePos :
		useLabel = useLabel + 'P'
	useLabel = useLabel + '_f'
	if useFeatPathZScore :
		useLabel = useLabel + 'Z'
	if useFeatTermWeights :
		useLabel = useLabel + 'T'
	if useFeatNeighbor :
		useLabel = useLabel + 'N'
	useLabel = useLabel + '_m{}'.format(negMultiplier)

	if printFlag :
		print("Using label: {}".format(useLabel))



	# 1) Load the gene-index dictionary & path names
	geneDict, pathDict = cl.getGeneAndPathDict(sDir)
	geneNames = list(geneDict.keys())
	geneNames.sort()
	pathNames = cl.removeInvertedPaths(pathDict)
	del pathDict



	# 2) Load the network general features
	numFN = 0
	if useFeatNeighbor :
		featNbVals, featNbNames = cl.getFeaturesNeighborhood(sDir, 'LogScale')
		featNbNames = np.ravel(featNbNames)
		numFN = len(featNbNames)
	#end if
	numTW = 0
	if useFeatTermWeights :
		featTWVals, featTWNames = cl.getFeaturesTerms(sDir, 'Orig')
		featTWNames = np.ravel(featTWNames)
		numTW = len(featTWNames)
	#end if



	# 3) Loop over the list of the sample subdirectories
	dSubDirs = cl.getSubDirectoryList(sDir)

	thisRound = 0
	for si in dSubDirs :
		thisRound += 1

		# Display directory to examine
		sv = si.split('/')
		if printFlag :
			print("\n{}/{}/".format(sv[-3],sv[-2]))

		# Create index lists for Known, Hidden, Unknown, TrueNeg from files
		giKnown, giUnknown, giHidden, giTrueNeg = cl.getGeneIndexLists(si, geneDict)
		giAll = list()
		giAll.extend(giKnown)
		giAll.extend(giUnknown)
		giAll.sort()


		# 4) Load the sample-specific features
		numFP = 0
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

		if useFeatPathZScore :
			if printFlag :
				print("    ... including path z-score features")
			features = np.hstack( (features, featZSVals) )
			featNames.extend(featZSNames)
		if useFeatNeighbor :
			if printFlag :
				print("    ... including neighborhood features")
			features = np.hstack( (features, featNbVals) )
			featNames.extend(np.ravel(featNbNames))
		if useFeatTermWeights :
			# Remove terms with no connection to gene set
			sumFTV = np.sum(featTWVals[giKnown,:], axis=0)
			keepIdx = np.nonzero(sumFTV)
			numTW = len(keepIdx[0])
			if printFlag :
				print("    ... including term membership features")
			features = np.hstack( (features, featTWVals[:,keepIdx[0]]) )
			featNames.extend(np.ravel(featTWNames[keepIdx]))
		# verify some features have been loaded
		if features.shape[1] == 0 :
			print("ERROR: No features were specified for classification.")
			sys.exit 
		#end if

		# Normalize the feature values
		features = cl.normalizeFeatureColumns(features)



		# Create the structure to rank the Unknown genes & paths
		geneRanks = np.zeros( (len(geneDict), 1), dtype=np.int32 )
		geneScores = np.zeros( (len(geneDict), 1), dtype=np.float32 )

#TODO: How to save feature rankings ??


		voteScores = np.zeros( (iterNumGenes, numVotes), dtype=np.float32)
#			voteScores = np.zeros( (len(geneDict), numVotes), dtype=np.float32)

		if printFlag :
			print("{} votes; cfier {}".format(numVotes, useCfier))
			print("  known: {}, total: {}, trainSet: {}".format( len(giKnown),
				iterNumGenes, (len(giKnown) * (1 + negMultiplier)) ))
		#end if



		# 6) Prepare the test/train vectors & labels
		# Extract the vectors for the pos sets

		retrySubSample = False
		retries = 0
		vote = 0
		while vote < numVotes :

			if retrySubSample :
				retrySubSample = False

				numSubSample = int(numSubSample * retrySubPortion) + 1
				retryIterKnown = random.sample(iterKnown, numSubSample)
				if len(retryIterKnown) < retryMinValid :
					retryIterKnown = random.sample(iterKnown, retryMinValid)

				posTrain = features[retryIterKnown,:]
				posTrainLabel = np.ones( (len(retryIterKnown), 1) ) * pLabel

				nExamples = min( negMultiplier * len(retryIterKnown), (iterNumGenes - len(retryIterKnown)))
			else :
				numSubSample = len(iterKnown)
	
				posTrain = features[iterKnown,:]
				posTrainLabel = np.ones( (len(iterKnown), 1) ) * pLabel

				nExamples = min( negMultiplier * len(iterKnown), (iterNumGenes - len(iterKnown)) )
			#end if

			# Extract the vectors for neg sets
			# as one-class: train with rand samp from Unknown
			#		test with all Unknown (TrueNeg + Hidden/TruePos)
			giTrainNeg = random.sample(iterUnknown, nExamples)
			negTrain = features[giTrainNeg,:]
			negTrainLabel = np.ones( (len(giTrainNeg), 1) ) * nLabel

			# Combine to create the full train & test data sets
			# as one-class:
			trainSet = np.vstack( (posTrain, negTrain) )
			trainLabel = np.vstack( (posTrainLabel, negTrainLabel) )

			testSet = features[giAll,:]

			# Some versions want the labels reshaped
			trainLabel = np.reshape(trainLabel, [trainLabel.shape[0],])


			# 7) Train classifier, predict on test, collect scores
			cfier = lm.LassoCV(alphas=useGivenRange, positive=usePos,
				max_iter=lMaxIter, normalize=lNorm, fit_intercept=lFitIcpt)
			cfier.fit(trainSet, trainLabel)
			foundAlpha = cfier.alpha_


			if printFlag :
				print("    Vote {} of {}; iters {:3d}, alpha {:.5f}, score {:.3f}; coeffs {}; sample {}".format(
					(vote + 1), numVotes, cfier.n_iter_, foundAlpha, cfier.score(trainSet, trainLabel),
					len(np.nonzero(cfier.coef_)[0]), len(posTrainLabel) ))
			#end if

			cfPredLabel = cfier.predict(testSet)
			cfPredLabel = np.ravel(cfPredLabel)

			# If no coeffs (train score == 0) try again
			if len(np.nonzero(cfier.coef_)[0]) <= 0 :
				if retries < (numVotes * 5) :
					retrySubSample = True
					vote = vote - 1
					retries += 1
				else :
					if printFlag :
						print("WARNING: used all retries.")
			else :
				numSubSample = len(iterKnown)
			#end if

			voteScores[:,vote] = cfPredLabel

			vote += 1
		#end loop (vote)


		# 8) Place the scores into the array and store across iterations

		# first, average across the random negative samples (votes)
#TODO: really, I should either normalize the score or vote across rank
		voteScores = mp.normalizeFeatureColumns(voteScores)
		voteAvgScore = np.mean(voteScores, axis=1)

		# then, place into full gene score array
		#	NOTE: carry the scores forward from each iteration
#			for u in range(len(iterUnknown)) :
#				geneScores[iterUnknown[u],iter] = voteAvgScore[u]
		for g in range(len(giAll)) :
			geneScores[giAll[g],1] = voteAvgScore[g]

#TODO: continue HERE


		# 10) Rank the genes across the iterations
#TODO: should I average these, or just take the last column ?
#	test that option later
		useScore = np.mean(geneScores[giUnknown,0:(itr+1)], axis=1)

		ranker = np.recarray(len(giUnknown),
			dtype=[('inverse', 'f4'), ('score', 'f4'), ('geneIdx', 'i4')])

		ranker['score'] = useScore
		ranker['inverse'] = np.multiply(useScore, -1)
		ranker['geneIdx'] = giUnknown

		ranker.sort(order=['inverse', 'geneIdx'])
	


		# 11) Output the ranked genes to file

		# write the file
		fname = 'ranked_genes-' + useLabel + '_Avg.txt'
		if printFlag :
			print("  Saving ranked genes to file {}".format(fname))
		with open(si+fname, 'w') as fout :
			firstRow = True
			for row in range(len(ranker)) :
				if not firstRow :
					fout.write('\n')
				fout.write('{:3.3f}{}{}'.format(ranker['score'][row],
					textDelim, geneNames[ranker['geneIdx'][row]]))
				firstRow = False
		#end with


		# 10-b) Rank the genes across the iterations
#TODO: should I average these, or just take the last column ?
#	test that option later
		useScore = geneScores[giUnknown,itr]

		ranker = np.recarray(len(giUnknown),
			dtype=[('inverse', 'f4'), ('score', 'f4'), ('geneIdx', 'i4')])
		ranker['score'] = useScore
		ranker['inverse'] = np.multiply(useScore, -1)
		ranker['geneIdx'] = giUnknown
		ranker.sort(order=['inverse', 'geneIdx'])
	
		# 11-b) Output the ranked genes to file
		# write the file
		fname = 'ranked_genes-' + useLabel + '_Last.txt'
		if printFlag :
			print("  Saving ranked genes to file {}".format(fname))
		with open(si+fname, 'w') as fout :
			firstRow = True
			for row in range(len(ranker)) :
				if not firstRow :
					fout.write('\n')
				fout.write('{:3.3f}{}{}'.format(ranker['score'][row],
					textDelim, geneNames[ranker['geneIdx'][row]]))
				firstRow = False
		#end with

		# 10-c) Rank the genes across the iterations
#TODO: should I average these, or just take the last column ?
#	test that option later
		useScore = geneScores[giUnknown,0]
		ranker = np.recarray(len(giUnknown),
			dtype=[('inverse', 'f4'), ('score', 'f4'), ('geneIdx', 'i4')])
		ranker['score'] = useScore
		ranker['inverse'] = np.multiply(useScore, -1)
		ranker['geneIdx'] = giUnknown
		ranker.sort(order=['inverse', 'geneIdx'])
			# 11-b) Output the ranked genes to file
		# write the file
		fname = 'ranked_genes-' + useLabel + '_First.txt'
		if printFlag :
			print("  Saving ranked genes to file {}".format(fname))
		with open(si+fname, 'w') as fout :
			firstRow = True
			for row in range(len(ranker)) :
				if not firstRow :
					fout.write('\n')
				fout.write('{:3.3f}{}{}'.format(ranker['score'][row],
					textDelim, geneNames[ranker['geneIdx'][row]]))
				firstRow = False
		#end with


		# 12) Output the selected feature info to file
#TODO: this


		# 13) Output the parameters to file
#TODO: this
		fname = 'parameters-' + useLabel + '.txt'
		with open(si+fname, 'w') as fout :
			fout.write('\n')
			fout.write('Sampling Method for Neg examples\n')
			fout.write('  as One-Class w/ iterations on the weaker predictions\n')
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
			elif useCfier == 2: 
				fout.write('method:{}ElasticNet\n'.format(textDelim))
			fout.write('positive:{}{}\n'.format(textDelim, usePos))
			# fout.write('alpha range:{}{}\n'.format(textDelim, useGivenRange))
			# fout.write('alpha chosen:{}{}\n'.format(textDelim, cfier.alpha_))
			fout.write('max_iter:{}{}\n'.format(textDelim, lMaxIter))
			fout.write('normalize:{}{}\n'.format(textDelim, lNorm))
			fout.write('fit_intercept:{}{}\n'.format(textDelim, lFitIcpt))
			fout.write('\n')
		#end with

		if printFlag :
			print("--{} of {}".format(thisRound, len(dSubDirs)))
	#end loop (si)

#end def ####### ####### ####### 




####### ####### ####### ####### 
# FUNCTION CALL

tstart = time.time()
print("\nPredicting gene ranks for {}".format(sDir))

# Main Function call
print("Calling function predictIterative() ...")
predictIterative(True)

print("--elapsed time: {:.3} (s)".format(time.time()-tstart))



print("\nDone.\n")