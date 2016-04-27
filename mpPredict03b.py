# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
#
# A different approach to ranking meta-paths.
# In this case, the collected counts between genes and the
#	training set will be normalized and used to weight
#	the metapaths by predictability, or importance. This
#	can be throught of as using feedback from the PathSim
#	metric to rank genes in order to influence and automate
#	the decision of which paths to use for prediction.
#
# This is the second part of the process: looking at the 
#	gene-gene similarity metric(s), apply LASSO regression
#	to select a small number of important paths, with
#	corresponding weights.
# ---------------------------------------------------------

import mpLibrary as mp
import time
import numpy as np
import gzip
from sklearn import linear_model as lm
import random



####### ####### ####### ####### 
# PARAMETERS

# Variables/Quantities
#percHide = [0, 10, 25, 33, 50]		# percent of genes to conceal


# Input names & locations
useNtwk = 1		# network & samples to use (0 means fake)

if useNtwk == 0 :
#	eName = 'fakeNtwk00_g2e3t10'
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
#	sPath = 'samplesFake/'
	dRoot = '../Dropbox/mp/outputFake/'
	dDir = 'pred03-batch-000'
else :
	eName = 'all_v1_g2e11t0'
	ePath = '../Dropbox/mp/networks/'
#	sPath = '../Dropbox/mp/samples-test1/'
	dRoot = '../Dropbox/mp/output/'
	dDir = 'pred03-batch-000'
#end if

# whether to write output files
writeOutput = True


# File names for similarity metrics
fSimilarity = 'SxySum.gz' 
	# SxySum.gz or Pxy.gz


# how to sample negative train/test set
sampleAsOneClass = False

# use .LassoCV vs .Lasso
useLCV = True

# LASSO params
if sampleAsOneClass :
	lAlpha = 0.005	# good for asOneClass
else :
	lAlpha = 0.0001 	# find good TwoClass alpha
#end if
lMaxIter = 100000	# .Lasso=1000, .LassoCV=100000
lNorm = True
lPos = False
lFitIcpt = True
lSelctn = 'random' # random vs cyclic

#TODO: will probably prefer using LassoCV to find alpha

# verbose feedback ?
newVerbose = True

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



# 2) Get a list of the sample subdirectories
dSubDirs = mp.getSubDirectoryList(dRoot+dDir)



# 3) For each sample (subdir), perform LASSO
#		save top paths & weights to a file

#for si in dSubDirs[4:6] :
for si in dSubDirs :

	# Display directory to examine
	sv = si.split('/')
	print("\n{}/{}/".format(sv[-3],sv[-2]))



	# ####### ####### ####### #######
	# 3a) Prepare the test/train vectors & labels
	print("  ... creating train & test data ...")

	# Read in the (modified) group PathSim matrix
	mp.setParamMatrixDT(np.float64)
	features = mp.readFileAsMatrix(si, fSimilarity)
	mp.setParamMatrixDT(-2)

	# Remove the mean, set L2 norm = 1
	featMean = np.mean(features, axis=0)
	features = np.subtract(features, featMean)
	featAbsMax = np.minimum(featMean, np.amax(features, axis=0))
	featAbsMax = np.add(featAbsMax, 1)	# hack so as not to / by 0
	features = np.divide(features, featAbsMax)

	# Create index lists for Known, Hidden, Unknown, TrueNeg
	gKnown = mp.readFileAsList(si+'known.txt')
	giKnown = mp.convertToIndices(gKnown, geneDict)
	gHidden = mp.readFileAsList(si+'concealed.txt')
	giHidden = mp.convertToIndices(gHidden, geneDict)
	giUnknown = [g for g in geneDict.values() if g not in giKnown]
	giTrueNeg = [g for g in giUnknown if g not in giHidden]
#	print len(gKnown), len(gHidden), len(giUnknown), len(giTrueNeg)

	# Extract the vectors for the pos sets
	posTrain = features[giKnown,:]
	posTest = features[giHidden,:]

	posTrainLabel = np.ones( (len(giKnown), 1) )
	posTestLabel = np.ones( (len(giHidden), 1) )


	# Extract the vectors for the neg sets
#TODO: best way to handle/produce negative train data ??
	if sampleAsOneClass == True :
		# as one-class: train with rand samp from Unknown
		#		test with all Unknown (TrueNeg + Hidden)
		nExamples = min( 2 * len(giKnown), (len(geneDict) - len(giKnown)) )
		giTrainNeg = random.sample(giUnknown, nExamples)

		negTrain = features[giTrainNeg]
		negTrainLabel = np.ones( (len(giTrainNeg), 1) ) * -1

		negTest = features[giTrueNeg]
		negTestLabel = np.ones( (len(giTrueNeg), 1) ) * -1

		giTest = np.hstack( (giHidden, giTrueNeg) )
	else :
		# as two-class: train rand samp from TrueNeg (1/2)
		#		test with remaining (TrueNeg * 1/2 + Hidden)
		nExamples = int( len(giTrueNeg) / 2 )
		giTrainNeg = random.sample(giTrueNeg, nExamples)
		giTestNeg = [g for g in giTrueNeg if g not in giTrainNeg]

		negTrain = features[giTrainNeg]
		negTrainLabel = np.ones( (len(giTrainNeg), 1) ) * -1

		negTest = features[giTestNeg]
		negTestLabel = np.ones( (len(giTestNeg), 1) ) * -1

		giTest = np.hstack( (giHidden, giTestNeg) )
	#end if


	# Combine to create the full train & test data setsa
	trainSet = np.vstack( (posTrain, negTrain) )
	trainLabel = np.vstack( (posTrainLabel, negTrainLabel) )

	testSet = np.vstack( (posTest, negTest) )
	testLabel = np.vstack( (posTestLabel, negTestLabel) )



	# ####### ####### ####### #######
	# 3b) Perform the regression analysis
	print("  ... performing regression ...")

	# Train LASSO	

	# Some versions want the labels reshaped
	trainLabel = np.reshape(trainLabel, [trainLabel.shape[0],])

	#TODO: exp w/ diff types: LassoCV, ElasticNet, MultiTaskElasticNet, MultiTaskLasso ..?
	if not useLCV :
		cfLasso = lm.Lasso(alpha=lAlpha, max_iter=lMaxIter, normalize=lNorm,
		 	positive=lPos, fit_intercept=lFitIcpt)#, selection=lSelctn)
		#TODO: Why did it start complaining about 'selection' ?
	else :
		cfLasso = lm.LassoCV(max_iter=lMaxIter, normalize=lNorm, positive=lPos,
			fit_intercept=lFitIcpt)
	#end if
	cfLasso.fit(trainSet, trainLabel)

	# The meaning of this score is questionable,
	#	mostly keeping it for curiosity
	cfScore = cfLasso.score(trainSet, trainLabel)
	print("On training data {}, score: {}".format(fSimilarity, cfScore))
	print("  using {} coefficients".format( len(np.nonzero(cfLasso.coef_)[0]) ))

	cfPredLabel = cfLasso.predict(testSet)



	# ####### ####### ####### #######
	# 3c) Output selected path results to file

	# Save the selected paths & scores/weights
	if writeOutput :

		textDelim = '\t'


		# Output the feature coefficients (mp weights)
		#	for the first metric
		cfCoefs = np.nonzero(cfLasso.coef_)[0]
		cfPaths = np.recarray( len(cfCoefs), dtype=[('path', 'i4'), ('weight', 'f4')] )
		row = 0
		for c in cfCoefs :
			cfPaths[row] = (c, cfLasso.coef_[c])
			row += 1
		cfPaths[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

		# write the file
		fPrefix = 'ranked_paths_Lasso-'+fSimilarity.rstrip('.txtgz')
		fname = mp.nameOutputFile(si, fPrefix)
		print("Saving top paths to file {}".format(fname))
		with open(si+fname, 'wb') as fout :
			fout.write('intercept:{}{}'.format(textDelim, cfLasso.intercept_))
#			fout.write('alpha:{0}{1}{0}max_iter:{0}{2}{0}'.format(textDelim, lAlpha, lMaxIter) +
#				'normalize:{0}{1}{0}positive:{0}{2}{0}'.format(textDelim, lNorm, lPos) +
#				'fit_intercept:{0}{1}{0}selection:{0}{2}{0}'.format(textDelim, lFitIcpt, lSelctn))
			for row in xrange(len(cfPaths)) :
				fout.write('\n{}{}{}'.format(cfPaths['weight'][row],
					textDelim, pathNames[cfPaths['path'][row]]))
		#end with
	#end if



	# ####### ####### ####### #######
	# 3d) Output ranked gene results to file

	# Save the genes from the test set to a file
	if writeOutput :

		textDelim = '\t'

#		geneList = geneDict.keys()
#		geneList.sort()

		# Sort the genes by (inverse) rank
		#	first metric
		cfGenes = np.recarray( len(cfPredLabel), dtype=[('gene', 'i4'), ('rank', 'f4')] )
		cfGenes['gene'] = giTest
		cfGenes['rank'] = cfPredLabel
		cfGenes[::-1].sort(order=['rank','gene'])

		# write the file
		fPrefix = 'ranked_genes-'+fSimilarity.rstrip('.txtgz')
		fname = mp.nameOutputFile(si, fPrefix)
		print("Saving ranked genes to file {}".format(fname))
		with open(si+fname, 'wb') as fout :
			firstRow = True
			for row in xrange(len(cfGenes)) :
				if not firstRow :
					fout.write('\n')
				fout.write('{:3.3f}{}{}'.format(cfGenes['rank'][row],
					textDelim, geneList[cfGenes['gene'][row]]))
				firstRow = False
		#end with
	#end if



	# ####### ####### ####### #######
	# 3e) Output parameters & scores to file

	# Save the parameters & results
	if writeOutput :

		textDelim = '\t'

		fname = mp.nameOutputFile(si, 'parameters')
		print("Saving params & stats to file {}".format(fname))
		with open(si+fname, 'wb') as fout :
			fout.write('\n')

			fout.write('Sampling Method for Neg examples\n')
			if sampleAsOneClass :
				fout.write('  as One-Class\n')
			else :
				fout.write('  as Two-Class\n')
			#end if
			fout.write('\n')

			fout.write('Lasso Parameters\n')
			if useLCV :
				fout.write('method:{}LassoCV (cross-validation)\n'.format(textDelim))
				fout.write('alpha:{}{}\n'.format(textDelim, cfLasso.alpha_))
			else :
				fout.write('method:{}Lasso (standard)\n'.format(textDelim))
				fout.write('alpha:{}{}\n'.format(textDelim, lAlpha))
			#end if
			fout.write('max_iter:{}{}\n'.format(textDelim, lMaxIter))
			fout.write('normalize:{}{}\n'.format(textDelim, lNorm))
			fout.write('positive:{}{}\n'.format(textDelim, lPos))
			fout.write('fit_intercept:{}{}\n'.format(textDelim, lFitIcpt))
			fout.write('selection:{}{}\n'.format(textDelim, lSelctn))
			fout.write('\n')

			fout.write('Similarity Metric:{}{}\n'.format(textDelim, fSimilarity))
			fout.write('Prediction Results\n')
			fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(np.nonzero(cfLasso.coef_)[0])))
			fout.write('Training score:{}{:3.3f}\n'.format(textDelim, cfLasso.score(trainSet, trainLabel)))
			fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, cfLasso.score(testSet, testLabel)))
			fout.write('\n')
		#end with
	#end if

#end loop



print("\nDone.\n")