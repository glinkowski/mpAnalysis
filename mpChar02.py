# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Characterization: define unique types of connections in set
#	Step 2:
#		find useful features & rank genes
#
# TODO: description
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
sDir = '../Dropbox/mp/output/char01-msigGood01'

# verbose feedback ?
verboseOutput = True

######## ######## ######## ######## 




######## ######## ######## ######## 
# PRIMARY FUNCTION
#TODO: add description
def predictIterative(printFlag) :

#TODO: remove / hardcode the following parameters
#	remove any superfluous ...

	# Parameters

	numVotes = 11
	usePos = True	
		# True/False: limit to only Positive coefficient values

	# Label for pos & neg labels
	pLabel = 1
	nLabel = 0
	negMultiplier = 1
	
	# LASSO params
	lMaxIter = 1000
	lNorm = True
	lFitIcpt = True

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

	# Control the iterations & error
	numVotes = 11
		# how many random samples for comparison
	retrySubPortion = 0.75
		# how many of Known to keep in new sub-sample
	retryMinValid = 9
		# minimum Known genes to use for PosTrain


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
		numFeatAll = len(featNames)
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


		voteScores = np.zeros( (len(giAll), numVotes), dtype=np.float32)
#			voteScores = np.zeros( (len(geneDict), numVotes), dtype=np.float32)

		if printFlag :
			print("{} votes; known: {}, total: {}, trainSet: {}".format( 
				numVotes, len(giKnown), len(giAll), (len(giKnown) * (1 + negMultiplier)) ))
		#end if



		# Store how many samples use certain features
		featT1List = np.zeros( (numFeatAll, 1), dtype=np.int16)
		featT5List = np.zeros( (numFeatAll, 1), dtype=np.int16)
		featTAList = np.zeros( (numFeatAll, 1), dtype=np.int16)
		# featT1Dict = dict()
		# featT1Set = set()
		# featT5Dict = dict()
		# featT5Set = set()
		# featTADict = dict()
		# featTASet = set()


		# 6) Prepare the test/train vectors & labels
		# Extract the vectors for the pos sets

		retrySubSample = False
		retries = 0
		vote = 0
		while vote < numVotes :

			if len(giKnown) < retryMinValid :
				retrySubSample = False

			if retrySubSample :
				retrySubSample = False

				numSubSample = int(numSubSample * retrySubPortion) + 1
				retryIterKnown = random.sample(giKnown, numSubSample)
				if len(retryIterKnown) < retryMinValid :
					retryIterKnown = random.sample(giKnown, retryMinValid)

				posTrain = features[retryIterKnown,:]
				posTrainLabel = np.ones( (len(retryIterKnown), 1) ) * pLabel

				nExamples = min( negMultiplier * len(retryIterKnown), (len(giAll) - len(retryIterKnown)))
			else :
				giKnown = giKnown
				numSubSample = len(giKnown)
	
				posTrain = features[giKnown,:]
				posTrainLabel = np.ones( (len(giKnown), 1) ) * pLabel

				nExamples = min( negMultiplier * len(giKnown), (len(giAll) - len(giKnown)) )
			#end if

			# Extract the vectors for neg sets
			# as one-class: train with rand samp from Unknown
			#		test with all Unknown (TrueNeg + Hidden/TruePos)
			giTrainNeg = random.sample(giUnknown, nExamples)
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

			# Else, collect info about the top-used features
			else :
				numSubSample = len(giKnown)

				# Extract indices corresponding to top 5 weighted features
				featWeights = cfier.coef_
				numFeats = len(np.nonzero(featWeights)[0])
				topFeats = np.ones( (numFeats), dtype=np.int32 ) * (-1)
				for num in range(numFeats) :
					featIdx = np.argmax(featWeights)
					topFeats[num] = featIdx
					featWeights[featIdx] = 0
				#end loop

				# Increment count for the Top 1 path
				featT1List[topFeats[0]] += 1

				# Increment count for the Top 5 paths
				for num in range(5) :
					if numFeats <= num :
						break
					#end if
					featT5List[topFeats[num]] += 1
				#end loop


				# Increment count for all non-zero paths
				for num in range(numFeats) :
					if numFeats <= num :
						break
					#end if
					featTAList[topFeats[num]] += 1
				#end loop
			#end if

			voteScores[:,vote] = cfPredLabel

			vote += 1
		#end loop (vote)


		# 8) Place the scores into the array and store across iterations

		# first, average across the random negative samples (votes)
#TODO: really, I should either normalize the score or vote across rank
		voteScores = cl.normalizeFeatureColumns(voteScores)
		voteAvgScore = np.mean(voteScores, axis=1)
		voteUnknownScore = voteAvgScore[giUnknown]

		ranker = np.recarray(len(giUnknown),
			dtype=[('inverse', 'f4'), ('score', 'f4'), ('geneIdx', 'i4')])

		ranker['score'] = voteUnknownScore
		ranker['inverse'] = np.multiply(voteUnknownScore, -1)
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
					'\t', geneNames[ranker['geneIdx'][row]]))
				firstRow = False
		#end with


		# 12) Output the selected feature info to file
#TODO: this

		# Sort the Top 1 paths
		featT1Sort = np.recarray( numFeatAll, dtype=[('featIdx', 'i4'), ('count', 'i4')])
		for row in range(numFeatAll) :
			featT1Sort['featIdx'][row] = row
			featT1Sort['count'][row] = featT1List[row]
		#end if
		featT1Sort[::-1].sort(order=['count', 'featIdx'])

		# Save the Top 1 paths to file
		fname = 'ranked_features_Top1-' + useLabel + '.txt'
		with open(si + fname, 'w') as fout :
			fout.write('Votes:{}{}'.format('\t', numVotes))
			row = 0
			nextVal = featT1Sort['count'][row]
			while nextVal != 0 :
				fout.write('\n{}{}{}'.format( nextVal, '\t',
					featNames[featT1Sort['featIdx'][row]]))
				row += 1
				nextVal = featT1Sort['count'][row]
		#end with

		# Sort the Top 5 paths
		featT5Sort = np.recarray( numFeatAll, dtype=[('featIdx', 'i4'), ('count', 'i4')])
		for row in range(numFeatAll) :
			featT5Sort['featIdx'][row] = row
			featT5Sort['count'][row] = featT5List[row]
		#end if
		featT5Sort[::-1].sort(order=['count', 'featIdx'])

		# Save the Top 5 paths to file
		fname = 'ranked_features_Top5-' + useLabel + '.txt'
		with open(si + fname, 'w') as fout :
			fout.write('Votes:{}{}'.format('\t', numVotes))
			row = 0
			nextVal = featT5Sort['count'][row]
			while nextVal != 0 :
				fout.write('\n{}{}{}'.format( nextVal, '\t',
					featNames[featT5Sort['featIdx'][row]]))
				row += 1
				nextVal = featT5Sort['count'][row]
		#end with

		# Sort the Top All Non-Zero paths
		featTASort = np.recarray( numFeatAll, dtype=[('featIdx', 'i4'), ('count', 'i4')])
		for row in range(numFeatAll) :
			featTASort['featIdx'][row] = row
			featTASort['count'][row] = featTAList[row]
		#end if
		featTASort[::-1].sort(order=['count', 'featIdx'])

		# Save the Top All Non-Zero paths to file
		fname = 'ranked_features_TopNZ-' + useLabel + '.txt'
		with open(si + fname, 'w') as fout :
			fout.write('Votes:{}{}'.format('\t', numVotes))
			row = 0
			nextVal = featTASort['count'][row]
			while nextVal != 0 :
				fout.write('\n{}{}{}'.format( nextVal, '\t',
					featNames[featTASort['featIdx'][row]]))
				row += 1
				nextVal = featTASort['count'][row]
		#end with



		# 13) Output the parameters to file
		fname = 'parameters-' + useLabel + '.txt'
		with open(si+fname, 'w') as fout :
			fout.write('\n')
			fout.write('Sampling Method for Neg examples\n')
			fout.write('  as One-Class w/ iterations on the weaker predictions\n')
			fout.write('\n')
			fout.write('Features Used\n')
			fout.write('path Z-Score:{}{}\n'.format('\t', useFeatPathZScore))
			fout.write('Neighborhood:{}{}\n'.format('\t', useFeatNeighbor))
			fout.write('Term Weights:{}{}\n'.format('\t', useFeatTermWeights))
			fout.write('\n')

#TODO: collect some stats (ie: common alphas, l1 ratios, etc)
			fout.write('Classifier Parameters\n')
			fout.write('method:{}Lasso\n'.format('\t'))
			fout.write('positive:{}{}\n'.format('\t', usePos))
			fout.write('alpha range:{}{}\n'.format('\t', useGivenRange))
			fout.write('alpha chosen:{}{}\n'.format('\t', cfier.alpha_))
			fout.write('max_iter:{}{}\n'.format('\t', lMaxIter))
			fout.write('normalize:{}{}\n'.format('\t', lNorm))
			fout.write('fit_intercept:{}{}\n'.format('\t', lFitIcpt))
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