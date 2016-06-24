# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Approach 2: learn top paths from PathSim sum, find genes
#	Version 2, Step 2(f)
#	use the 2-class version of sampling
#	lasso, lasso positive, elastic net positive
#
# The first approach tries to first find the most
#	important metapaths, then rank genes by similarity
#	along those. This approach instead finds the set
#	similarity along each metapath, then tries to learn
#	the most important paths.
#
# Step 2f: From the pre-built per-gene feature vector,
#	apply LASSO linear regression to select a small number
#	of important paths. Follow these paths to find genes
#	connected to the set, and rank those by similarity.
# NOTE: This version uses the two-class method of defining
#	negatives, where it is assumed the true neg and true pos
#	are fully defined and known beforehand. Then
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

retCutoffs = [50, 100, 200, 500, 1000, 2000]


# LASSO params
lAlpha02 = [0.0008, 0.0007, 0.0006, 0.0005, 0.0004, 0.0003]
lAlphaPos = [0.001, 0.0007, 0.0006, 0.0005, 0.0004, 0.0003, 0.0001, 0.00008]
lMaxIter = 10000
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
print("\nPerforming regression(s) on {} ...".format(sDir))

mp.setParamVerbose(newVerbose)



# 1) Load the gene-index dictionary & path names
geneDict, pathDict = mp.getGeneAndPathDict(sDir)
geneNames = list(geneDict.keys())
geneNames.sort()
pathNames = mp.removeInvertedPaths(pathDict)
del pathDict


# 2) Loop over all the subdirectories
dSubDirs = mp.getSubDirectoryList(sDir)

thisRound = 0
#for si in dSubDirs[0:1] :
for si in dSubDirs :

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

#TODO: add neighborhood features ?

	# Normalize & Create train/test sets
	trainSet, trainLabel, testSet, testLabel, giTest = mp.createTrainTestSets(si, geneDict, features, False)



	# ####### ####### ####### #######
	# 4) Perform the regression analysis, Lasso 2-class
	print("Lasso Regression ...")

	# Train Lasso Regression
#	cfLassoPos = lm.LassoCV(alphas=lAlpha02, max_iter=lMaxIter, normalize=lNorm,
#		positive=False, fit_intercept=lFitIcpt)	
	cfLassoPos = lm.LassoCV(max_iter=lMaxIter, normalize=lNorm,
		positive=False, fit_intercept=lFitIcpt)
	cfLassoPos.fit(trainSet, trainLabel)

	# The meaning of this score is questionable,
	#	mostly keeping it for curiosity
	cfScore = cfLassoPos.score(trainSet, trainLabel)
	print("Lasso Regression score: {}".format(cfScore))
	print("  coefficients: {}".format( len(np.nonzero(cfLassoPos.coef_)[0]) ))
	print("  chosen alpha: {}".format( cfLassoPos.alpha_ ))
	print("  iterations: {}".format( cfLassoPos.n_iter_ ))
#	print("  using {} coefficients, alpha = {}".format(
#		len(np.nonzero(cfLassoPos.coef_)[0]), cfLassoPos.alpha_ ))

	cfPredLabel = cfLassoPos.predict(testSet)
	cfPredLabel = np.ravel(cfPredLabel)


	# 5) Output results to file, Lasso 2-class

	# Save the selected paths & scores/weights
	# 	feature coefficients are the metapath weights
	cfCoefs = np.nonzero(cfLassoPos.coef_)[0]
	cfPaths = np.recarray( len(cfCoefs), dtype=[('path', 'i4'), ('weight', 'f4')] )
	for row in range(len(cfCoefs)) :
		cfPaths[row] = (row, cfCoefs[row])
	cfPaths[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

	# write the file
	fname = 'ranked_paths-Lasso_2ClassStd.txt'
	print("Saving data for the Lasso standard approach ...")
	print("  Saving top paths to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		fout.write('intercept:{}{}'.format(textDelim, cfLassoPos.intercept_))
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
		fout.write('alpha:{}{}\n'.format(textDelim, cfLassoPos.alpha_))
		fout.write('max_iter:{}{}\n'.format(textDelim, lMaxIter))
		fout.write('normalize:{}{}\n'.format(textDelim, lNorm))
		fout.write('positive:{}{}\n'.format(textDelim, 'False'))
		fout.write('fit_intercept:{}{}\n'.format(textDelim, lFitIcpt))
		fout.write('\n')

		fout.write('Similarity Metric:{}PathSim sum over set\n'.format(textDelim))
		fout.write('Prediction Results\n')
		fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(np.nonzero(cfLassoPos.coef_)[0])))
		fout.write('Training score:{}{:3.3f}\n'.format(textDelim, cfLassoPos.score(trainSet, trainLabel)))
		fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, cfLassoPos.score(testSet, testLabel)))
		fout.write('\n')
	#end with


	# ####### ####### ####### #######
	# 4) Perform the regression analysis, Lasso 2-class
	print("Lasso (positive) Regression ...")

	# Train Lasso Regression
#	cfLassoPos = lm.LassoCV(alphas=lAlphaPos, max_iter=lMaxIter, normalize=lNorm,
#		positive=False, fit_intercept=lFitIcpt)
	cfLassoPos = lm.LassoCV(max_iter=lMaxIter, normalize=lNorm,
		positive=True, fit_intercept=lFitIcpt)
	cfLassoPos.fit(trainSet, trainLabel)

	# The meaning of this score is questionable,
	#	mostly keeping it for curiosity
	cfScore = cfLassoPos.score(trainSet, trainLabel)
	print("Lasso Regression score: {}".format(cfScore))
#	print("  using {} coefficients, alpha = {}".format(
#		len(np.nonzero(cfLassoPos.coef_)[0]), cfLassoPos.alpha_ ))
	print("  coefficients: {}".format( len(np.nonzero(cfLassoPos.coef_)[0]) ))
	print("  chosen alpha: {}".format( cfLassoPos.alpha_ ))
	print("  iterations: {}".format( cfLassoPos.n_iter_ ))

	cfPredLabel = cfLassoPos.predict(testSet)
	cfPredLabel = np.ravel(cfPredLabel)


	# 5) Output results to file, Lasso 2-class

	# Save the selected paths & scores/weights
	# 	feature coefficients are the metapath weights
	cfCoefs = np.nonzero(cfLassoPos.coef_)[0]
	cfPaths = np.recarray( len(cfCoefs), dtype=[('path', 'i4'), ('weight', 'f4')] )
	for row in range(len(cfCoefs)) :
		cfPaths[row] = (row, cfCoefs[row])
	cfPaths[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

	# write the file
	fname = 'ranked_paths-Lasso_2C_Pos.txt'
	print("Saving data for the Lasso standard approach ...")
	print("  Saving top paths to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		fout.write('intercept:{}{}'.format(textDelim, cfLassoPos.intercept_))
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
	fname = 'ranked_genes-Lasso_2C_Pos.txt'
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
	fname = 'parameters-Lasso_2C_Pos.txt'
	with open(si+fname, 'w') as fout :
		fout.write('\n')
		fout.write('Sampling Method for Neg examples\n')
		fout.write('  as Two-Class\n')
		fout.write('\n')

		fout.write('Lasso Parameters\n')
		fout.write('method:{}Lasso (pos coeffs)\n'.format(textDelim))
		fout.write('a_range:{}{}\n'.format(textDelim, lAlpha02))
		fout.write('alpha:{}{}\n'.format(textDelim, cfLassoPos.alpha_))
		fout.write('max_iter:{}{}\n'.format(textDelim, lMaxIter))
		fout.write('normalize:{}{}\n'.format(textDelim, lNorm))
		fout.write('positive:{}{}\n'.format(textDelim, 'True'))
		fout.write('fit_intercept:{}{}\n'.format(textDelim, lFitIcpt))
		fout.write('\n')

		fout.write('Similarity Metric:{}PathSim sum over set\n'.format(textDelim))
		fout.write('Prediction Results\n')
		fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(np.nonzero(cfLassoPos.coef_)[0])))
		fout.write('Training score:{}{:3.3f}\n'.format(textDelim, cfLassoPos.score(trainSet, trainLabel)))
		fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, cfLassoPos.score(testSet, testLabel)))
		fout.write('\n')
	#end with



	# ####### ####### ####### #######
	# 4) Perform the regression analysis, Elastic Net (Pos Coeffs only)
	print("Elastic Net w/ only Positive coefficients ...")

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

	cfPredLabel = cfENet.predict(testSet)
	cfPredLabel = np.ravel(cfPredLabel)



	# 5) Output results to file, Elastic Net (Pos Coeffs only)

	# Save the selected paths & scores/weights
	# 	feature coefficients are the metapath weights
	cfCoefs = cfENet.coef_
	cfPaths = np.recarray( len(cfCoefs), dtype=[('path', 'i4'), ('weight', 'f4')] )
	for row in range(len(cfCoefs)) :
		cfPaths[row] = (row, cfCoefs[row])
	cfPaths[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

	# write the file
	fname = 'ranked_paths-ElasticNet_2C_Pos.txt'
	print("Saving data for the Elastic Net approach ...")
	print("  Saving top paths to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		fout.write('intercept:{}{}'.format(textDelim, cfENet.intercept_))
		for row in range(len(cfPaths)) :
			fout.write('\n{}{}{}'.format(cfPaths['weight'][row],
				textDelim, pathNames[cfPaths['path'][row]]))
	#end with

	#Sort the genes by (inverse) rank
	cfGenes = np.recarray( len(cfPredLabel), dtype=[('gene', 'i4'), ('rank', 'f4')] )
	cfGenes['gene'] = giTest
	cfGenes['rank'] = np.ravel(cfPredLabel)
	cfGenes[::-1].sort(order=['rank','gene'])

	# write the file
	fname = 'ranked_genes-ElasticNet_2C_Pos.txt'
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
	fname = 'parameters-ElasticNet_2C_Pos.txt'
	with open(si+fname, 'w') as fout :
		fout.write('\n')
		fout.write('Sampling Method for Neg examples\n')
		fout.write('  as Two-Class\n')
		fout.write('\n')

		fout.write('Log Regression Parameters\n')
		fout.write('method:{}ElasticNet CV\n'.format(textDelim))
		fout.write('ratio range:{}{}\n'.format(textDelim, enNAlphas))
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