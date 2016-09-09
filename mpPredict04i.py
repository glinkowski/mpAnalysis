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


# verbose feedback ?
verbose = True



#TODO: remove / hardcode the following parameters


# Adjustible classifier parameters
useCfier = 1
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
useFeatTermWeights = True
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


# text delimiter in output files
textDelim = '\t'

####### ####### ####### ####### 

#TODO: If score == 0,
# find a new alpha (findAlpha flag)
# select new random set


####### ####### ####### ####### 
# ANCILLARY FUNCTIONS

def getGeneIndexLists(path, gDict) :

	# Create index lists for Known, Hidden, Unknown, TrueNeg
	gKnown = mp.readFileAsList(path + 'known.txt')
	giKnown = mp.convertToIndices(gKnown, gDict)
	gHidden = mp.readFileAsList(path + 'concealed.txt')
	giHidden = mp.convertToIndices(gHidden, gDict)
	giUnknown = [g for g in gDict.values() if g not in giKnown]
	giTrueNeg = [g for g in giUnknown if g not in giHidden]

	return giKnown, giUnknown, giHidden, giTrueNeg
#end def ####### ####### ####### 


####### ####### ####### ####### 




####### ####### ####### ####### 
# PRIMARY FUNCTION

def predictIterative(printFlag) :

	if printFlag :
		print("\nPerforming regression(s) on {}".format(sDir))

	mp.setParamVerbose(verbose)


	# 0) Create the useLabel variable
	# string: label for the output files
	# ie: ClusVote_c<Las/Enet/Log/SVM><P for Pos>_f<P for pathsim><Z for z-score>
	#	<T for term weights><N for neighborhood>
	useLabel = 'Iter{}V{}_c'.format(numIterations, numVotes)
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
	useLabel = useLabel + '_m{}'.format(negMultiplier)
	if retryOnZeroCoeffs :
		useLabel = useLabel + '_wRS' # indicating resample on 0 score
	#end if

	if printFlag :
		print("Using label: {}".format(useLabel))


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
	if printFlag :
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
	dSubDirs = mp.getSubDirectoryList(sDir)

	thisRound = 0
	#for si in dSubDirs[0:1] :
	for si in dSubDirs :
		thisRound += 1

		# Display directory to examine
		sv = si.split('/')
		if printFlag :
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
			if printFlag :
				print("    ... including PathSim sum features")
			features = np.hstack( (features, featPSVals) )
			featNames.extend(featPSNames)
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
		features = mp.normalizeFeatureColumns(features)



		# Peform N recursive iterations, each voted across K random samples
		
		# Create the structure to rank the Unknown genes & paths
		geneRanks = np.zeros( (len(geneDict), 1), dtype=np.int32 )
		geneScores = np.zeros( (len(geneDict), numIterations), dtype=np.float32 )

#TODO: How to save feature rankings ??

		# set the gene indices for the first iteration
		iterKnown = giKnown
		iterUnknown = giUnknown
		iterAll = list()
		iterAll.extend(iterKnown)
		iterAll.extend(iterUnknown)
		iterAll.sort()

		for itr in range(numIterations) :

			# store the results for each random sample
#			iterNumGenes = len(iterKnown) + len(iterUnknown)
			iterNumGenes = len(iterAll)
			voteScores = np.zeros( (iterNumGenes, numVotes), dtype=np.float32)
#			voteScores = np.zeros( (len(geneDict), numVotes), dtype=np.float32)

			if printFlag :
				print("  iteration {} of {}; {} votes; cfier {}".format( (itr + 1),
					numIterations, numVotes, useCfier))
				print("  known: {}, total: {}, trainSet: {}".format( len(iterKnown),
					iterNumGenes, (len(iterKnown) * (1 + negMultiplier)) ))
			#end if


			# 6) Prepare the test/train vectors & labels
			# Extract the vectors for the pos sets
			posTrain = features[iterKnown,:]
			posTrainLabel = np.ones( (len(iterKnown), 1) ) * pLabel

			findNewAlpha = True
			retries = 0
			for vote in range(numVotes) :

				# Extract the vectors for neg sets
				# as one-class: train with rand samp from Unknown
				#		test with all Unknown (TrueNeg + Hidden/TruePos)
				nExamples = min( negMultiplier * len(iterKnown), (iterNumGenes - len(iterKnown)) )
				giTrainNeg = random.sample(iterUnknown, nExamples)
				negTrain = features[giTrainNeg,:]
				negTrainLabel = np.ones( (len(giTrainNeg), 1) ) * nLabel

				# Combine to create the full train & test data sets
				# as one-class:
				trainSet = np.vstack( (posTrain, negTrain) )
				trainLabel = np.vstack( (posTrainLabel, negTrainLabel) )

#				testSet = features[iterUnknown,:]
				testSet = features[iterAll,:]

				# Some versions want the labels reshaped
				trainLabel = np.reshape(trainLabel, [trainLabel.shape[0],])


				# 7) Train classifier, predict on test, collect scores

#TODO: add other classifier options ??
				if useCfier == 1 :	# 1 = Lasso
					if findNewAlpha :
						cfier = lm.LassoCV(alphas=useGivenRange, positive=usePos,
							max_iter=lMaxIter, normalize=lNorm, fit_intercept=lFitIcpt)
						cfier.fit(trainSet, trainLabel)
						foundAlpha = cfier.alpha_
						findNewAlpha = False
					else :
						cfier = lm.Lasso(alpha=foundAlpha, max_iter=lMaxIter, normalize=lNorm,
							positive=usePos, fit_intercept=lFitIcpt)
						cfier.fit(trainSet, trainLabel)
					#end if
				else :
					print("ERROR: specified classifier unrecognized: useCfier = {}".format(useCfier))
				#end if


#TODO: Decide on verbose vs printFlag
				if verbose :
					# view quick statistics from this training session
					if useCfier < 3 :
						if printFlag :
							print("    Vote {}-{}; iters {:3d}, alpha {:.5f}, score {:.3f}; coeffs {}".format((itr + 1), (vote + 1),
								cfier.n_iter_, foundAlpha, cfier.score(trainSet, trainLabel), len(np.nonzero(cfier.coef_)[0])))
							if useCfier == 2 :	# 2 = ElasticNet
								print("    l1 ratio: {}".format( cfier.l1_ratio_ ))
				#end if

#TODO: What if 0 coefficients are chosen? Select new sample? Let it slide (affect avg?)
				if useCfier < 3 :
					cfPredLabel = cfier.predict(testSet)
				#end if
				cfPredLabel = np.ravel(cfPredLabel)

				# If no coeffs (train score == 0) try again
				if retryOnZeroCoeffs :
					if len(np.nonzero(cfier.coef_)[0]) <= 0 :
						if retries < (numVotes * 3) :
							findNewAlpha = True
							vote -= 1
							retries += 1
				#end if

				voteScores[:,vote] = cfPredLabel
			#end loop (vote)

			# In case an iteration is no longer useful (all scores == 0)
			if (vote > 0) and (cfier.score(trainSet, trainLabel) <= 0) :
				break


			# 8) Place the scores into the array and store across iterations

			# first, average across the random negative samples (votes)
#TODO: really, I should either normalize the score or vote across rank
#	Does the value of the intersection matter here?
			voteScores = mp.normalizeFeatureColumns(voteScores)
			voteAvgScore = np.mean(voteScores, axis=1)

			# then, place into full gene score array
			#	NOTE: carry the scores forward from each iteration
#			for u in range(len(iterUnknown)) :
#				geneScores[iterUnknown[u],iter] = voteAvgScore[u]
			for g in range(len(iterAll)) :
				geneScores[iterAll[g],itr] = voteAvgScore[g]
			for i in range(itr + 1, numIterations) :
				geneScores[:,i] = geneScores[:,itr]
			#end loop


			# 9) Select Known & Unknown for the next round
#TODO: Base this on mis-labelled Known genes
#	for now, just take a percentage of least-confident scores

			# find the cutoff value for scores to keep
#			idxKeep = len(iterAll) - int(len(iterAll) / float(numIterations))
			cutoffIdx = iterNumGenes - int(iterNumGenes / float(numIterations))
			absScore = np.absolute(voteAvgScore)
			absScore.sort()
			cutoffVal = absScore[cutoffIdx]

			# get the upper & lower indices to extract for next round
			# x = 0
			# testVal = voteAvgScore
			# while testVal >= cutoffVal :

			# extract indices for any genes scoring less than cutoff
			iterKeep = list()
			for x in range(len(iterAll)) :
				if abs(voteAvgScore[x]) < cutoffVal :
					iterKeep.append(iterAll[x])
			#end loop

			# find intersections of Keep w/ previous Known & Unknown
			setKeep = set(iterKeep)
			newKnown = [gi for gi in iterKnown if gi in setKeep]
#			newKnown.sort()
			newUnknown = [gi for gi in iterUnknown if gi in setKeep]
#			newUnknown.sort()

			# set the gene indices for the next iteration
			iterKnown = newKnown
			iterUnknown = newUnknown
			iterAll = list()
			iterAll.extend(iterKnown)
			iterAll.extend(iterUnknown)
			iterAll.sort()

			numKnown = len(iterKnown)
			numUnknown = len(iterUnknown)
			if (numKnown <= numExitKnown) or (numUnknown <= numExitUnknown) :
				if printFlag :
					print("known: {}, unknown: {}; exiting loop".format(numKnown, numUnknown))
				break
		#end loop (itr)

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