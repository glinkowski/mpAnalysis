# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Approach 2: learn top paths from PathSim sum, find genes
#	Version 2, Step 2(e)
#	log regression; elastic net
#
# The first approach tries to first find the most
#	important metapaths, then rank genes by similarity
#	along those. This approach instead finds the set
#	similarity along each metapath, then tries to learn
#	the most important paths.
#
# Step 2b: From the pre-built per-gene feature vector,
#	apply LASSO linear regression to select a small number
#	of important paths. Follow these paths to find genes
#	connected to the set, and rank those by similarity.
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
dDir = 'pred04-set01'
dRoot = '../Dropbox/mp/output/'

# File name containing feature vectors
fSimilarity = 'features_PathSim.gz' 

retCutoffs = [50, 100, 200, 500, 1000, 2000]


# Log Regression params
#lgCs = np.logspace(-4, 1, 6)
lgCs = 3
#lgCs = 10
lgPenalty = 'l2'
lgDual = False
lgMaxIter = 1000

# Elastic Net params
#enRatios = [0.3, 0.5, 0.8, 0.95]
enRatios = [0.2, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99]
#enRatios = [0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999]
enNAlphas = 11
#enNAlphas = 25
enMaxIter = 100000
#enPos = False
enFitIncept = True
enNorm = True
enCopy = True


# Lables for pos & neg training sets
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
if dDir.endswith('/') :
	dPath = dRoot + dDir
else :
	dPath = dRoot + dDir + '/'
with open(dPath + 'parameters.txt', 'r') as fin :
	line = fin.readline()

	line = fin.readline()
	line = line.rstrip()
	lv = line.split(textDelim)
	eName = lv[1]

	line = fin.readline()
	line = line.rstrip()
	lv = line.split(textDelim)
	ePath = lv[1]
#end with

print("Creating the gene-index dictionary.")
geneDict = mp.readGenesFile(ePath, eName)
geneList = list(geneDict.keys())
geneList.sort()

print("Reading in the path names.")
pathDict = mp.readKeyFile(ePath, eName)
pathNames = mp.removeInvertedPaths(pathDict)
del pathDict



# 2) Loop over all the subdirectories
dSubDirs = mp.getSubDirectoryList(dRoot+dDir)

thisRound = 0
for si in dSubDirs[0:1] :
#for si in dSubDirs :

	# Display directory to examine
	sv = si.split('/')
	print("\n{}/{}/".format(sv[-3],sv[-2]))



	# ####### ####### ####### #######
	# 3) Prepare the test/train vectors & labels
	print("  ... creating train & test data ...")

	# Read in the feature vectors
	features = mp.readFileAsMatrix(si, fSimilarity)
	# NOTE: previous version of mpPredict04a had extra zeros at end of vectors
	#	discarding those columns
	features = features[:,0:len(pathNames)]


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
#	print("{}, {}, {}, {}".format(len(gKnown), len(gHidden), len(giUnknown), len(giTrueNeg)))


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



	# ####### ####### ####### #######
	# 4) Perform the regression analysis, Logarithmic
	print("Logarithmic Regression ...")
	print("  ... performing regression ...")

	# Train Log Regression
	cfLog = lm.LogisticRegressionCV(Cs=lgCs, penalty=lgPenalty, dual=lgDual, max_iter=lgMaxIter)
	cfLog.fit(trainSet, trainLabel)

	# The meaning of this score is questionable,
	#	mostly keeping it for curiosity
	cfScore = cfLog.score(trainSet, trainLabel)
	print("Log Regression score: {}".format(cfScore))
	print("  using {} coefficients, C = {} = e^{}".format(
		len(np.nonzero(cfLog.coef_)[0]), cfLog.C_,  np.log(cfLog.C_) ))

	cfPredLabel = cfLog.predict(testSet)
#	cfPredLabel = cfPredLabel.reshape( (len(cfPredLabel), 1) )
	cfPredLabel = np.ravel(cfPredLabel)



	# 5) Output results to file, Logarithmic

	# Save the selected paths & scores/weights
	# 	feature coefficients are the metapath weights
	cfCoefs = cfLog.coef_[0]
#	print(cfLog.coef_)
#	print(cfCoefs)
	cfPaths = np.recarray( len(cfCoefs), dtype=[('path', 'i4'), ('weight', 'f4')] )
	for row in range(len(cfCoefs)) :
		cfPaths[row] = (row, cfCoefs[row])
	# row = 0
	# for c in cfCoefs :
	# 	cfPaths[row] = (c, cfLasso01.coef_[c])
	# 	row += 1
	cfPaths[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

	# write the file
	fname = 'ranked_paths-LogRegression.txt'
	print("Saving data for the Log Regression approach ...")
	print("  Saving top paths to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		fout.write('intercept:{}{}'.format(textDelim, cfLog.intercept_))
		for row in range(len(cfPaths)) :
			fout.write('\n{}{}{}'.format(cfPaths['weight'][row],
				textDelim, pathNames[cfPaths['path'][row]]))
	#end with

#	# Save the genes from the test set to a file
#	mp.writeRankedGenes02(si, 'LogRegression', cfPredLabel, geneDict, giKnown,
#		giHidden, retCutoffs)

	#Sort the genes by (inverse) rank
	cfGenes = np.recarray( len(cfPredLabel), dtype=[('gene', 'i4'), ('rank', 'f4')] )
	cfGenes['gene'] = giTest
#	cfGenes['rank'] = np.ravel(cfPredLabel)
	cfGenes['rank'] = cfPredLabel
	cfGenes[::-1].sort(order=['rank','gene'])

	# write the file
	fname = 'ranked_genes-LogRegression.txt'
	print("  Saving ranked genes to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		firstRow = True
		for row in range(len(cfGenes)) :
			if not firstRow :
				fout.write('\n')
			fout.write('{:3.3f}{}{}'.format(cfGenes['rank'][row],
				textDelim, geneList[cfGenes['gene'][row]]))
			firstRow = False
	#end with

	# Save the parameters & results
	fname = 'parameters-LogRegression.txt'
	with open(si+fname, 'w') as fout :
		fout.write('\n')
		fout.write('Sampling Method for Neg examples\n')
		fout.write('  as One-Class\n')
		fout.write('\n')

		fout.write('Log Regression Parameters\n')
		fout.write('method:{}Log Reg CV\n'.format(textDelim))
		fout.write('C range:{}{}\n'.format(textDelim, lgCs))
		fout.write('C chosen:{}{}\n'.format(textDelim, cfLog.C_))
		fout.write('n_iter:{}{}\n'.format(textDelim, cfLog.n_iter_))
		fout.write('intercept:{}{}\n'.format(textDelim, cfLog.intercept_))
		fout.write('\n')

		fout.write('Similarity Metric:{}PathSim sum over set\n'.format(textDelim))
		fout.write('Prediction Results\n')
		fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(np.nonzero(cfLog.coef_)[0])))
		fout.write('Training score:{}{:3.3f}\n'.format(textDelim, cfLog.score(trainSet, trainLabel)))
		fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, cfLog.score(testSet, testLabel)))
		fout.write('\n')
	#end with



	# ####### ####### ####### #######
	# 4) Perform the regression analysis, Elastic Net (non-Pos)
	print("Elastic Net ...")
	print("  ... performing regression ...")

	# Train Log Regression
	cfENet = lm.ElasticNetCV(l1_ratio=enRatios, positive=False, fit_intercept=enFitIncept,
		n_alphas=enNAlphas, normalize=enNorm, copy_X=enCopy, max_iter=enMaxIter)
	cfENet.fit(trainSet, trainLabel)

	# The meaning of this score is questionable,
	#	mostly keeping it for curiosity
	cfScore = cfENet.score(trainSet, trainLabel)
	print("Elastic Net stuff:".format(cfScore))
	print("  # coefficients: {}".format( len(np.nonzero(cfENet.coef_)[0]) ))
	print("  l1 ratio: {}".format(cfENet.l1_ratio_))
	print("  alpha: {}".format(cfENet.alpha_))
	print("  iterations: {}".format(cfENet.n_iter_))
#	print("  all alphas: {}".format(cfENet.alphas_))

	cfPredLabel = cfENet.predict(testSet)
#	cfPredLabel = cfPredLabel.reshape( (len(cfPredLabel), 1) )
	cfPredLabel = np.ravel(cfPredLabel)



	# 5) Output results to file, Elastic Net (non-Pos)

	# Save the selected paths & scores/weights
	# 	feature coefficients are the metapath weights
	cfCoefs = cfENet.coef_
#	print(cfENet.coef_)
#	print(cfCoefs)
	cfPaths = np.recarray( len(cfCoefs), dtype=[('path', 'i4'), ('weight', 'f4')] )
	for row in range(len(cfCoefs)) :
		cfPaths[row] = (row, cfCoefs[row])
	# row = 0
	# for c in cfCoefs :
	# 	cfPaths[row] = (c, cfLasso01.coef_[c])
	# 	row += 1
	cfPaths[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

	# write the file
	fname = 'ranked_paths-ElasticNet.txt'
	print("Saving data for the Elastic Net approach ...")
	print("  Saving top paths to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		fout.write('intercept:{}{}'.format(textDelim, cfENet.intercept_))
		for row in range(len(cfPaths)) :
			fout.write('\n{}{}{}'.format(cfPaths['weight'][row],
				textDelim, pathNames[cfPaths['path'][row]]))
	#end with

#	# Save the genes from the test set to a file
#	mp.writeRankedGenes02(si, 'LogRegression', cfPredLabel, geneDict, giKnown,
#		giHidden, retCutoffs)

	#Sort the genes by (inverse) rank
	cfGenes = np.recarray( len(cfPredLabel), dtype=[('gene', 'i4'), ('rank', 'f4')] )
	cfGenes['gene'] = giTest
	cfGenes['rank'] = np.ravel(cfPredLabel)
	cfGenes[::-1].sort(order=['rank','gene'])

	# write the file
	fname = 'ranked_genes-ElasticNet.txt'
	print("  Saving ranked genes to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		firstRow = True
		for row in range(len(cfGenes)) :
			if not firstRow :
				fout.write('\n')
			fout.write('{:3.3f}{}{}'.format(cfGenes['rank'][row],
				textDelim, geneList[cfGenes['gene'][row]]))
			firstRow = False
	#end with

	# Save the parameters & results
	fname = 'parameters-ElasticNet.txt'
	with open(si+fname, 'w') as fout :
		fout.write('\n')
		fout.write('Sampling Method for Neg examples\n')
		fout.write('  as One-Class\n')
		fout.write('\n')

		fout.write('Log Regression Parameters\n')
		fout.write('method:{}ElasticNet CV\n'.format(textDelim))
		fout.write('ratio range:{}{}\n'.format(textDelim, lgCs))
		fout.write('ratio chosen:{}{}\n'.format(textDelim, cfENet.l1_ratio_))
		fout.write('alpha chosen:{}{}\n'.format(textDelim, cfENet.alpha_))
		fout.write('n_iter:{}{}\n'.format(textDelim, cfENet.n_iter_))
		fout.write('intercept:{}{}\n'.format(textDelim, cfENet.intercept_))
		fout.write('\n')

		fout.write('Similarity Metric:{}PathSim sum over set\n'.format(textDelim))
		fout.write('Prediction Results\n')
		fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(np.nonzero(cfENet.coef_)[0])))
		fout.write('Training score:{}{:3.3f}\n'.format(textDelim, cfENet.score(trainSet, trainLabel)))
		fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, cfENet.score(testSet, testLabel)))
		fout.write('\n')
	#end with



	# ####### ####### ####### #######
	# 4) Perform the regression analysis, Elastic Net (Pos Coeffs only)
	print("Elastic Net w/ only Positive coefficients ...")
	print("  ... performing regression ...")

	# Train Log Regression
	cfENet = lm.ElasticNetCV(l1_ratio=enRatios, positive=True, fit_intercept=enFitIncept,
		n_alphas=enNAlphas, normalize=enNorm, copy_X=enCopy, max_iter=enMaxIter)
	cfENet.fit(trainSet, trainLabel)

	# The meaning of this score is questionable,
	#	mostly keeping it for curiosity
	cfScore = cfENet.score(trainSet, trainLabel)
	print("Elastic Net stuff:".format(cfScore))
	print("  # coefficients: {}".format( len(np.nonzero(cfENet.coef_)[0]) ))
	print("  l1 ratio: {}".format(cfENet.l1_ratio_))
	print("  alpha: {}".format(cfENet.alpha_))
	print("  iterations: {}".format(cfENet.n_iter_))
#	print("  all alphas: {}".format(cfENet.alphas_))

	cfPredLabel = cfENet.predict(testSet)
#	cfPredLabel = cfPredLabel.reshape( (len(cfPredLabel), 1) )
	cfPredLabel = np.ravel(cfPredLabel)



	# 5) Output results to file, Elastic Net (Pos Coeffs only)

	# Save the selected paths & scores/weights
	# 	feature coefficients are the metapath weights
	cfCoefs = cfENet.coef_
#	print(cfENet.coef_)
#	print(cfCoefs)
	cfPaths = np.recarray( len(cfCoefs), dtype=[('path', 'i4'), ('weight', 'f4')] )
	for row in range(len(cfCoefs)) :
		cfPaths[row] = (row, cfCoefs[row])
	# row = 0
	# for c in cfCoefs :
	# 	cfPaths[row] = (c, cfLasso01.coef_[c])
	# 	row += 1
	cfPaths[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

	# write the file
	fname = 'ranked_paths-ElasticNet_Pos.txt'
	print("Saving data for the Elastic Net approach ...")
	print("  Saving top paths to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		fout.write('intercept:{}{}'.format(textDelim, cfENet.intercept_))
		for row in range(len(cfPaths)) :
			fout.write('\n{}{}{}'.format(cfPaths['weight'][row],
				textDelim, pathNames[cfPaths['path'][row]]))
	#end with

#	# Save the genes from the test set to a file
#	mp.writeRankedGenes02(si, 'LogRegression', cfPredLabel, geneDict, giKnown,
#		giHidden, retCutoffs)

	#Sort the genes by (inverse) rank
	cfGenes = np.recarray( len(cfPredLabel), dtype=[('gene', 'i4'), ('rank', 'f4')] )
	cfGenes['gene'] = giTest
	cfGenes['rank'] = np.ravel(cfPredLabel)
	cfGenes[::-1].sort(order=['rank','gene'])

	# write the file
	fname = 'ranked_genes-ElasticNet_Pos.txt'
	print("  Saving ranked genes to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		firstRow = True
		for row in range(len(cfGenes)) :
			if not firstRow :
				fout.write('\n')
			fout.write('{:3.3f}{}{}'.format(cfGenes['rank'][row],
				textDelim, geneList[cfGenes['gene'][row]]))
			firstRow = False
	#end with

	# Save the parameters & results
	fname = 'parameters-ElasticNet_Pos.txt'
	with open(si+fname, 'w') as fout :
		fout.write('\n')
		fout.write('Sampling Method for Neg examples\n')
		fout.write('  as One-Class\n')
		fout.write('\n')

		fout.write('Log Regression Parameters\n')
		fout.write('method:{}ElasticNet CV\n'.format(textDelim))
		fout.write('ratio range:{}{}\n'.format(textDelim, lgCs))
		fout.write('ratio chosen:{}{}\n'.format(textDelim, cfENet.l1_ratio_))
		fout.write('alpha chosen:{}{}\n'.format(textDelim, cfENet.alpha_))
		fout.write('n_iter:{}{}\n'.format(textDelim, cfENet.n_iter_))
		fout.write('intercept:{}{}\n'.format(textDelim, cfENet.intercept_))
		fout.write('\n')

		fout.write('Similarity Metric:{}PathSim sum over set\n'.format(textDelim))
		fout.write('Prediction Results\n')
		fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(np.nonzero(cfENet.coef_)[0])))
		fout.write('Training score:{}{:3.3f}\n'.format(textDelim, cfENet.score(trainSet, trainLabel)))
		fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, cfENet.score(testSet, testLabel)))
		fout.write('\n')
	#end with



	thisRound += 1
	print("    --{} of {}".format(thisRound, len(dSubDirs)))
	print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))
#end loop



#TODO: append output file in root of batch ?



print("\nDone.\n")