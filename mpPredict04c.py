# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Approach 2: learn top paths from PathSim sum, find genes
#	Version 2, Step 2(c)
#	LASSO + Voting
#
# The first approach tries to first find the most
#	important metapaths, then rank genes by similarity
#	along those. This approach instead finds the set
#	similarity along each metapath, then tries to learn
#	the most important paths.
#
# Step 2c: From the pre-built per-gene feature vector,
#	apply LASSO linear regression to select a small number
#	of important paths. Follow these paths to find genes
#	connected to the set, and rank those by similarity.
# NEW: Do this several times, and vote among the gene
#	rankings to create a final rank
# ---------------------------------------------------------

import mpLibrary as mp
import time
import numpy as np
import gzip
from sklearn import linear_model as lm
import random



####### ####### ####### ####### 
# PARAMETERS

# number of rounds over which to vote
numVotes = 100

# folder containing the pre-processed samples
dDir = 'pred04-set01'

# Input names & locations
useNtwk = 1		# network & samples to use (0 means fake)
if useNtwk == 0 :
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
	dRoot = 'outputFake/'
else :
	eName = 'all_v3beta_g2e9t0'
	ePath = '../Dropbox/mp/networks/'
	dRoot = '../Dropbox/mp/output/'
#end if

# File name containing feature vectors
fSimilarity = 'features_PathSim.gz' 


# LASSO params
lAlphas = [0.005, 0.003, 0.001, 0.0007, 0.0004, 0.0001]
lMaxIter = 10000
lNorm = True
lPos = False
lFitIcpt = True
#lSelctn = 'random' # random vs cyclic

pLabel = 1
nLabel = 0


# verbose feedback ?
newVerbose = True

textDelim = '\t'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()
print("")

mp.setParamVerbose(newVerbose)



# 1) Load the gene-index dictionary & path names
print("Creating the gene-index dictionary.")
geneDict = mp.readGenesFile(ePath, eName)
geneList = list(geneDict.keys())
geneList.sort()

print("Reading in the path names.")
pathDict = mp.readKeyFile(ePath, eName)
pathNames = mp.removeInvertedPaths(pathDict)
del pathDict



# 2) Loop over the list of the sample subdirectories
dSubDirs = mp.getSubDirectoryList(dRoot+dDir)

thisRound = 0
#for si in dSubDirs[0:1] :
for si in dSubDirs :

	# Display directory to examine
	sv = si.split('/')
	print("\n{}/{}/".format(sv[-3],sv[-2]))

#	subSubDir = sv[-2]
#	sv = subSubDir.split('-')
#	print("  hidden: {}%, sample: {}".format(sv[0], sv[1]))


	# 3) Create the structure to rank the Unknown genes & paths

	# Create index lists for Known, Hidden, Unknown, TrueNeg
	gKnown = mp.readFileAsList(si+'known.txt')
	giKnown = mp.convertToIndices(gKnown, geneDict)
	gHidden = mp.readFileAsList(si+'concealed.txt')
	giHidden = mp.convertToIndices(gHidden, geneDict)
	giUnknown = [g for g in geneDict.values() if g not in giKnown]
#	print(giUnknown[0:21])
#	giUnknown.sort()	# just want to be certain of the order
	giTrueNeg = [g for g in giUnknown if g not in giHidden]
#	print("{}, {}, {}, {}".format(len(gKnown), len(gHidden), len(giUnknown), len(giTrueNeg)))

	# tradeoffs: dict vs array ... loop now or loop later?
#	gRankSums = dict()
#	for gi in giUnknown :
#		gRankSums[gi] = 0
	gRankSums = np.zeros( len(geneList) )
	pValSums = np.zeros( len(pathNames) )


	# ####### ####### ####### #######
	# 4) Use LassoCV to find best alpha for this sample & fold

	# 4a) Prepare the test/train vectors & labels
	print("  ... creating train & test data ...")

	# Read in the feature vectors
	features = mp.readFileAsMatrix(si, fSimilarity)

	# NOTE: previous version of mpPredict04a had extra zeros at end of vectors
	#	discarding those columns
	features = features[:,0:len(pathNames)]

	# Remove the mean; set L2 norm = 1
	featMean = np.mean(features, axis=0)
	features = np.subtract(features, featMean)
	featAbsMax = np.minimum(featMean, np.amax(features, axis=0))
	featAbsMax = np.add(featAbsMax, 1)	# hack so as not to / by 0
	features = np.divide(features, featAbsMax)

	# Extract the vectors for the pos sets
	posTrain = features[giKnown,:]
	posTest = features[giHidden,:]

	posTrainLabel = np.ones( (len(giKnown), 1) ) * pLabel
	posTestLabel = np.ones( (len(giHidden), 1) ) * pLabel

	# Extract the vectors for neg & Test sets
	# as one-class: train with rand samp from Unknown
	#		test with all Unknown (TrueNeg + Hidden)
	nExamples = min( 2 * len(giKnown), (len(geneDict) - len(giKnown)) )
	giTrainNeg = random.sample(giUnknown, nExamples)

	negTrain = features[giTrainNeg]
	negTrainLabel = np.ones( (len(giTrainNeg), 1) ) * nLabel

	negTest = features[giTrueNeg]
	negTestLabel = np.ones( (len(giTrueNeg), 1) ) * nLabel

	giTest = np.hstack( (giHidden, giTrueNeg) )

	# Combine to create the full train & test data sets
	# as one-class:
	trainSet = np.vstack( (posTrain, negTrain) )
	trainLabel = np.vstack( (posTrainLabel, negTrainLabel) )

	testSet = np.vstack( (posTest, negTest) )
	testLabel = np.vstack( (posTestLabel, negTestLabel) )

	# Some versions want the labels reshaped
	trainLabel = np.reshape(trainLabel, [trainLabel.shape[0],])
	testLabel = np.reshape(testLabel, [testLabel.shape[0],])


	# 4b) Perform the regression analysis
	print("  ... performing first regression ...")

	# Train LASSO, as one-class	
	cfLCV = lm.LassoCV(alphas=lAlphas, max_iter=lMaxIter, normalize=lNorm,
		positive=lPos, fit_intercept=lFitIcpt)

	cfLCV.fit(trainSet, trainLabel)

	# save the alpha for use in later regressions
	useAlpha = cfLCV.alpha_

	# The meaning of this score is questionable,
	#	mostly keeping it for curiosity
	cfScore = cfLCV.score(trainSet, trainLabel)
	print("One-Class score: {}".format(cfScore))
	print("  using {} coefficients, alpha = {}".format(
		len(np.nonzero(cfLCV.coef_)[0]), useAlpha ))

	cfPredLabel = cfLCV.predict(testSet)

	# Save the metapath (coefficient) values
#	print(features.shape)
#	print("{}, {}".format(len(pValSums), len(cfLCV.coef_)))
	pValSums = np.add(pValSums, cfLCV.coef_)

	# Sort the genes by (inverse) score
	cfGenes = np.recarray( len(cfPredLabel), dtype=[('gene', 'i4'), ('score', 'f4')] )
	cfGenes['gene'] = giTest	#giTest is same as giUnknown for one-class approach
	cfGenes['score'] = cfPredLabel
	cfGenes[::-1].sort(order=['score','gene'])

	# Add rank to gRankSums
	rank = 0
	for gi in cfGenes['gene'] :
		rank += 1
		gRankSums[gi] += rank
	#end loop

	# keep track of the average intercept
	avgIncept = cfLCV.intercept_
	avgNumCoeffs = len(np.nonzero(cfLCV.coef_)[0])
	avgTrainScore = cfScore
	avgTestScore = cfLCV.score(testSet, testLabel)



	# ####### ####### ####### #######
	# 5) Repeat Lasso process N times using best alpha
	print(" ... performing {} regressions ...".format( (numVotes - 1) ))
	for n in range(numVotes - 1) : 

		# 5a) Build train/test sets (pos data doesn't change)

		# Extract the vectors for neg & Test sets
		# as one-class: train with rand samp from Unknown
		#		test with all Unknown (TrueNeg + Hidden)
		nExamples = min( 2 * len(giKnown), (len(geneDict) - len(giKnown)) )
		giTrainNeg = random.sample(giUnknown, nExamples)

		negTrain = features[giTrainNeg]
		negTrainLabel = np.ones( (len(giTrainNeg), 1) ) * nLabel

		negTest = features[giTrueNeg]
		negTestLabel = np.ones( (len(giTrueNeg), 1) ) * nLabel

		giTest = np.hstack( (giHidden, giTrueNeg) )

		# Combine to create the full train & test data sets
		# as one-class:
		trainSet = np.vstack( (posTrain, negTrain) )
		trainLabel = np.vstack( (posTrainLabel, negTrainLabel) )

		testSet = np.vstack( (posTest, negTest) )
		testLabel = np.vstack( (posTestLabel, negTestLabel) )

		# Some versions want the labels reshaped
		trainLabel = np.reshape(trainLabel, [trainLabel.shape[0],])
		testLabel = np.reshape(testLabel, [testLabel.shape[0],])


		# 5b) Perform the regression analysis
		# Train LASSO, as one-class	
		cfLPlain= lm.Lasso(alpha=useAlpha, max_iter=lMaxIter, normalize=lNorm,
			positive=lPos, fit_intercept=lFitIcpt)

		cfLPlain.fit(trainSet, trainLabel)
		# Save the metapath (coefficient) values
		pValSums = np.add(pValSums, cfLPlain.coef_)

		# The meaning of this score is questionable,
		#	mostly keeping it for curiosity
		cfScore = cfLPlain.score(trainSet, trainLabel)
		if newVerbose :
			print("    score: {:1.3f}, # ceoffs: {}".format(cfScore, len(np.nonzero(cfLPlain.coef_)[0])))
#		elif not (thisRound % 25) :
#			print("    completed {} of {}".format(thisRound, (numVotes - 1) ))
	
		cfPredLabel = cfLPlain.predict(testSet)

		# Sort the genes by (inverse) score
		cfGenes = np.recarray( len(cfPredLabel), dtype=[('gene', 'i4'), ('score', 'f4')] )
		cfGenes['gene'] = giTest	#giTest is same as giUnknown for one-class approach
		cfGenes['score'] = cfPredLabel
		cfGenes[::-1].sort(order=['score','gene'])

		# Add rank to gRankSums
		rank = 0
		for gi in cfGenes['gene'] :
			rank += 1
			gRankSums[gi] += rank
		#end loop

		# keep track of some stats
		avgIncept += cfLPlain.intercept_
		avgNumCoeffs += len(np.nonzero(cfLPlain.coef_)[0])
		avgTrainScore += cfScore
		avgTestScore += cfLPlain.score(testSet, testLabel)

	#end loop



	# ####### ####### ####### #######
	# 6) Output results to file(s)

	# Save the selected paths & scores/weights (that are non-zero)
	# 	feature coefficients are the metapath weights
	cfPaths = np.recarray( len(np.nonzero(pValSums)[0]), dtype=[('path', 'i4'), ('weight', 'f4')] )
	row = 0
	for pi in range(len(pathNames)) :
		if pValSums[pi] != 0.0 :
			cfPaths[row] = (pi, pValSums[pi])
			row += 1
	cfPaths[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

	# write the file
	fname = 'ranked_paths-Lasso_Voting1C_{}.txt'.format(numVotes)
	print("Saving data for the Voting (One-Class) approach ...")
	print("  Saving top paths to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		fout.write('avg intercept:{}{}'.format( textDelim, (avgIncept / float(numVotes)) ))
		for row in range(len(cfPaths)) :
			fout.write('\n{}{}{}'.format(cfPaths['weight'][row],
				textDelim, pathNames[cfPaths['path'][row]]))
	#end with


	# Save the ranked genes to a file
	# Sort the genes by rank sum
	cfGenes = np.recarray( len(giUnknown), dtype=[('gene', 'i4'), ('rank', 'f4')] )
	row = 0
	for gi in giUnknown :
		cfGenes[row] = (gi, gRankSums[gi])
		row += 1
#	cfGenes['gene'] = giUnknown
#	cfGenes['rank'] = gRankSums
	cfGenes.sort(order=['rank','gene'])

	# write the file
	fname = 'ranked_genes-Lasso_Voting1C_{}.txt'.format(numVotes)
	print("  Saving ranked genes to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		firstRow = True
		for row in range(len(cfGenes)) :
			if not firstRow :
				fout.write('\n')
			fout.write('{:3.3f}{}{}'.format( (cfGenes['rank'][row] / float(numVotes)),
				textDelim, geneList[cfGenes['gene'][row]] ))
			firstRow = False
	#end with


	# Save the parameters & results
	fname = 'parameters-Lasso_Voting1C_{}.txt'.format(numVotes)
	with open(si+fname, 'w') as fout :
		fout.write('\n')
		fout.write('Sampling Method for Neg examples\n')
		fout.write('  as One-Class\n')
		fout.write('\n')

		fout.write('Lasso Parameters\n')
		fout.write('method:{}Lasso (standard)\n'.format(textDelim))
		fout.write('alpha range:{}{}\n'.format(textDelim, lAlphas))
		fout.write('alpha chosen:{}{}\n'.format(textDelim, useAlpha))
		fout.write('max_iter:{}{}\n'.format(textDelim, lMaxIter))
		fout.write('normalize:{}{}\n'.format(textDelim, lNorm))
		fout.write('positive:{}{}\n'.format(textDelim, lPos))
		fout.write('fit_intercept:{}{}\n'.format(textDelim, lFitIcpt))
		fout.write('\n')

		fout.write('Similarity Metric:{}PathSim sum over set\n'.format(textDelim))
		fout.write('Prediction Results\n')
		fout.write('nonzero coefficients:{}{}\n'.format(textDelim, (avgNumCoeffs / float(numVotes)) ))
		fout.write('Training score:{}{:3.3f}\n'.format(textDelim, (avgTrainScore / float(numVotes)) ))
		fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, (avgTestScore / float(numVotes)) ))
		fout.write('\n')
	#end with

	thisRound += 1
	print("    --{} of {}".format(thisRound, len(dSubDirs)))
	print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))
#end loop



print("\nDone.\n")