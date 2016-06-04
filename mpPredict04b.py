# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Approach 2: learn top paths from PathSim sum, find genes
#	Version 2, Step 2(b)
#	LASSO linear regression
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

# Input names & locations
useNtwk = 1		# network & samples to use (0 means fake)
if useNtwk == 0 :
#	eName = 'fakeNtwk01_g3e4t1'
#	ePath = 'networks/'
	dRoot = 'outputFake/'
else :
#	eName = 'all_v3beta_g2e9t0'
#	ePath = '../Dropbox/mp/networks/'
	dRoot = '../Dropbox/mp/output/'
#end if

# File name containing feature vectors
fSimilarity = 'features_PathSim.gz' 


# LASSO params
lAlpha01 = 0.0006
lAlpha02 = 0.00002
lAlphaPos = 0.005
lMaxIter = 10000
lNorm = True
lPos = False
lFitIcpt = True
lSelctn = 'random' # random vs cyclic

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



# 2) Get a list of the sample subdirectories
dSubDirs = mp.getSubDirectoryList(dRoot+dDir)



# 3) For each sample (subdir), perform LASSO
#		save top paths & weights to a file

#for si in dSubDirs[0:1] :
for si in dSubDirs :

	# Display directory to examine
	sv = si.split('/')
	print("\n{}/{}/".format(sv[-3],sv[-2]))



	# ####### ####### ####### #######
	# 3a) Prepare the test/train vectors & labels
	print("  ... creating train & test data ...")

	# Read in the feature vectors
#	mp.setParamMatrixDT(np.float64)
	features = mp.readFileAsMatrix(si, fSimilarity)
#	mp.setParamMatrixDT(-2)


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
	nExamples01 = min( 2 * len(giKnown), (len(geneDict) - len(giKnown)) )
	giTrainNeg01 = random.sample(giUnknown, nExamples01)

	negTrain01 = features[giTrainNeg01]
	negTrainLabel01 = np.ones( (len(giTrainNeg01), 1) ) * nLabel

	negTest01 = features[giTrueNeg]
	negTestLabel01 = np.ones( (len(giTrueNeg), 1) ) * nLabel

	giTest01 = np.hstack( (giHidden, giTrueNeg) )


	# Extract the vectors for neg & Test sets
	# as two-class: train rand samp from TrueNeg (1/2)
	#		test with remaining (TrueNeg * 1/2 + Hidden)
	nExamples02 = int( len(giTrueNeg) / 2 )
	giTrainNeg02 = random.sample(giTrueNeg, nExamples02)
	giTestNeg02 = [g for g in giTrueNeg if g not in giTrainNeg02]

	negTrain02 = features[giTrainNeg02]
	negTrainLabel02 = np.ones( (len(giTrainNeg02), 1) ) * nLabel

	negTest02 = features[giTestNeg02]
	negTestLabel02 = np.ones( (len(giTestNeg02), 1) ) * nLabel

	giTest02 = np.hstack( (giHidden, giTestNeg02) )


	# Combine to create the full train & test data sets
	# as one-class:
	trainSet01 = np.vstack( (posTrain, negTrain01) )
	trainLabel01 = np.vstack( (posTrainLabel, negTrainLabel01) )

	testSet01 = np.vstack( (posTest, negTest01) )
	testLabel01 = np.vstack( (posTestLabel, negTestLabel01) )

	# Some versions want the labels reshaped
	trainLabel01 = np.reshape(trainLabel01, [trainLabel01.shape[0],])
	testLabel01 = np.reshape(testLabel01, [testLabel01.shape[0],])


	# as two-class:
	trainSet02 = np.vstack( (posTrain, negTrain02) )
	trainLabel02 = np.vstack( (posTrainLabel, negTrainLabel02) )

	testSet02 = np.vstack( (posTest, negTest02) )
	testLabel02 = np.vstack( (posTestLabel, negTestLabel02) )

	# Some versions want the labels reshaped
	trainLabel02 = np.reshape(trainLabel02, [trainLabel02.shape[0],])
	testLabel02 = np.reshape(testLabel02, [testLabel02.shape[0],])



	# ####### ####### ####### #######
	# 3b) Perform the regression analysis
	print("  ... performing regression ...")

	# Train LASSO, as one-class	
	cfLasso01 = lm.Lasso(alpha=lAlpha01, max_iter=lMaxIter, normalize=lNorm,
		 	positive=lPos, fit_intercept=lFitIcpt)#, selection=lSelctn)
	#end if
	cfLasso01.fit(trainSet01, trainLabel01)

	# The meaning of this score is questionable,
	#	mostly keeping it for curiosity
	cfScore01 = cfLasso01.score(trainSet01, trainLabel01)
	print("One-Class score: {}".format(cfScore01))
	print("  using {} coefficients".format( len(np.nonzero(cfLasso01.coef_)[0]) ))

	cfPredLabel01 = cfLasso01.predict(testSet01)


	# ####### ####### ####### #######
	# 3c) Output results to file, One-Class

	# Save the selected paths & scores/weights
	# 	feature coefficients are the metapath weights
	cfCoefs = np.nonzero(cfLasso01.coef_)[0]
	cfPaths = np.recarray( len(cfCoefs), dtype=[('path', 'i4'), ('weight', 'f4')] )
	row = 0
	for c in cfCoefs :
		cfPaths[row] = (c, cfLasso01.coef_[c])
		row += 1
	cfPaths[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

	# write the file
#	fPrefix = 'ranked_paths-Lasso_1Class'
#	fname = mp.nameOutputFile(si, fPrefix)
	fname = 'ranked_paths-Lasso_1Class.txt'
	print("Saving data for the One-Class approach ...")
	print("  Saving top paths to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		fout.write('intercept:{}{}'.format(textDelim, cfLasso01.intercept_))
		for row in range(len(cfPaths)) :
			fout.write('\n{}{}{}'.format(cfPaths['weight'][row],
				textDelim, pathNames[cfPaths['path'][row]]))
	#end with


	# Save the genes from the test set to a file
	# Sort the genes by (inverse) rank
	cfGenes = np.recarray( len(cfPredLabel01), dtype=[('gene', 'i4'), ('rank', 'f4')] )
	cfGenes['gene'] = giTest01
	cfGenes['rank'] = cfPredLabel01
	cfGenes[::-1].sort(order=['rank','gene'])

	# write the file
#	fPrefix = 'ranked_genes-Lasso_1Class'
#	fname = mp.nameOutputFile(si, fPrefix)
	fname = 'ranked_genes-Lasso_1Class.txt'
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
#	fname = mp.nameOutputFile(si, 'parameters-Lasso_1Class')
#	print("  Saving params & stats to file {}".format(fname))
	fname = 'parameters-Lasso_1Class.txt'
	with open(si+fname, 'w') as fout :
		fout.write('\n')
		fout.write('Sampling Method for Neg examples\n')
		fout.write('  as One-Class\n')
		fout.write('\n')

		fout.write('Lasso Parameters\n')
		fout.write('method:{}Lasso (standard)\n'.format(textDelim))
		fout.write('alpha:{}{}\n'.format(textDelim, lAlpha01))
		fout.write('max_iter:{}{}\n'.format(textDelim, lMaxIter))
		fout.write('normalize:{}{}\n'.format(textDelim, lNorm))
		fout.write('positive:{}{}\n'.format(textDelim, lPos))
		fout.write('fit_intercept:{}{}\n'.format(textDelim, lFitIcpt))
		fout.write('selection:{}{}\n'.format(textDelim, lSelctn))
		fout.write('\n')

		fout.write('Similarity Metric:{}PathSim sum over set\n'.format(textDelim))
		fout.write('Prediction Results\n')
		fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(np.nonzero(cfLasso01.coef_)[0])))
		fout.write('Training score:{}{:3.3f}\n'.format(textDelim, cfLasso01.score(trainSet01, trainLabel01)))
		fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, cfLasso01.score(testSet01, testLabel01)))
		fout.write('\n')
	#end with


	# ####### ####### ####### #######
	# 3b) Perform the regression analysis

	# Train LASSO, as two-class	
	cfLasso02 = lm.Lasso(alpha=lAlpha02, max_iter=lMaxIter, normalize=lNorm,
		 	positive=lPos, fit_intercept=lFitIcpt)#, selection=lSelctn)
	#end if
	cfLasso02.fit(trainSet02, trainLabel02)

	# The meaning of this score is questionable,
	#	mostly keeping it for curiosity
	cfScore02 = cfLasso02.score(trainSet02, trainLabel02)
	print("Two-Class score: {}".format(cfScore02))
	print("  using {} coefficients".format( len(np.nonzero(cfLasso02.coef_)[0]) ))

	cfPredLabel02 = cfLasso02.predict(testSet02)


	# 3c) Output results to file, Two-Class

	# Save the selected paths & scores/weights
	# 	feature coefficients are the metapath weights
	cfCoefs = np.nonzero(cfLasso02.coef_)[0]
	cfPaths = np.recarray( len(cfCoefs), dtype=[('path', 'i4'), ('weight', 'f4')] )
	row = 0
	for c in cfCoefs :
		cfPaths[row] = (c, cfLasso02.coef_[c])
		row += 1
	cfPaths[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

	# write the file
#	fPrefix = 'ranked_paths-Lasso_2Class'
#	fname = mp.nameOutputFile(si, fPrefix)
	fname = 'ranked_paths-Lasso_2Class.txt'
	print("Saving data for the Two-Class approach ...")
	print("  Saving top paths to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		fout.write('intercept:{}{}'.format(textDelim, cfLasso02.intercept_))
		for row in range(len(cfPaths)) :
			fout.write('\n{}{}{}'.format(cfPaths['weight'][row],
				textDelim, pathNames[cfPaths['path'][row]]))
	#end with


	# Save the genes from the test set to a file
	# Sort the genes by (inverse) rank
	cfGenes = np.recarray( len(cfPredLabel02), dtype=[('gene', 'i4'), ('rank', 'f4')] )
	cfGenes['gene'] = giTest02
	cfGenes['rank'] = cfPredLabel02
	cfGenes[::-1].sort(order=['rank','gene'])

	# write the file
#	fPrefix = 'ranked_genes-Lasso_2Class'
#	fname = mp.nameOutputFile(si, fPrefix)
	fname = 'ranked_genes-Lasso_2Class.txt'
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
#	fname = mp.nameOutputFile(si, 'parameters-Lasso_2Class')
#	print("  Saving params & stats to file {}".format(fname))
	fname = 'parameters-Lasso_2Class.txt'
	with open(si+fname, 'w') as fout :
		fout.write('\n')
		fout.write('Sampling Method for Neg examples\n')
		fout.write('  as Two-Class\n')
		fout.write('\n')

		fout.write('Lasso Parameters\n')
		fout.write('method:{}Lasso (standard)\n'.format(textDelim))
		fout.write('alpha:{}{}\n'.format(textDelim, lAlpha02))
		fout.write('max_iter:{}{}\n'.format(textDelim, lMaxIter))
		fout.write('normalize:{}{}\n'.format(textDelim, lNorm))
		fout.write('positive:{}{}\n'.format(textDelim, lPos))
		fout.write('fit_intercept:{}{}\n'.format(textDelim, lFitIcpt))
		fout.write('selection:{}{}\n'.format(textDelim, lSelctn))
		fout.write('\n')

		fout.write('Similarity Metric:{}PathSim sum over set\n'.format(textDelim))
		fout.write('Prediction Results\n')
		fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(np.nonzero(cfLasso02.coef_)[0])))
		fout.write('Training score:{}{:3.3f}\n'.format(textDelim, cfLasso02.score(trainSet02, trainLabel02)))
		fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, cfLasso02.score(testSet02, testLabel02)))
		fout.write('\n')
	#end with



	# ####### ####### ####### #######
	# 3b) Perform the regression analysis, Positive coeffs

	# Train LASSO, as one-class	
	cfLasso01 = lm.Lasso(alpha=lAlphaPos, max_iter=lMaxIter, normalize=lNorm,
		 	positive=True, fit_intercept=lFitIcpt)#, selection=lSelctn)
	#end if
	cfLasso01.fit(trainSet01, trainLabel01)

	# The meaning of this score is questionable,
	#	mostly keeping it for curiosity
	cfScore01 = cfLasso01.score(trainSet01, trainLabel01)
	print("Pos-Coeffs score: {}".format(cfScore01))
	print("  using {} coefficients".format( len(np.nonzero(cfLasso01.coef_)[0]) ))

	cfPredLabel01 = cfLasso01.predict(testSet01)


	# ####### ####### ####### #######
	# 3c) Output results to file, One-Class

	# Save the selected paths & scores/weights
	# 	feature coefficients are the metapath weights
	cfCoefs = np.nonzero(cfLasso01.coef_)[0]
	cfPaths = np.recarray( len(cfCoefs), dtype=[('path', 'i4'), ('weight', 'f4')] )
	row = 0
	for c in cfCoefs :
		cfPaths[row] = (c, cfLasso01.coef_[c])
		row += 1
	cfPaths[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

	# write the file
#	fPrefix = 'ranked_paths-Lasso_Pos'
#	fname = mp.nameOutputFile(si, fPrefix)
	fname = 'ranked_paths-Lasso_Pos.txt'
	print("Saving data for the Positive approach ...")
	print("  Saving top paths to file {}".format(fname))
	with open(si+fname, 'w') as fout :
		fout.write('intercept:{}{}'.format(textDelim, cfLasso01.intercept_))
		for row in range(len(cfPaths)) :
			fout.write('\n{}{}{}'.format(cfPaths['weight'][row],
				textDelim, pathNames[cfPaths['path'][row]]))
	#end with


	# Save the genes from the test set to a file
	# Sort the genes by (inverse) rank
	cfGenes = np.recarray( len(cfPredLabel01), dtype=[('gene', 'i4'), ('rank', 'f4')] )
	cfGenes['gene'] = giTest01
	cfGenes['rank'] = cfPredLabel01
	cfGenes[::-1].sort(order=['rank','gene'])

	# write the file
#	fPrefix = 'ranked_genes-Lasso_Pos'
#	fname = mp.nameOutputFile(si, fPrefix)
	fname = 'ranked_genes-Lasso_Pos.txt'
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
#	fname = mp.nameOutputFile(si, 'parameters-Lasso_Pos')
#	print("  Saving params & stats to file {}".format(fname))
	fname = 'parameters-Lasso_Pos.txt'
	with open(si+fname, 'w') as fout :
		fout.write('\n')
		fout.write('Sampling Method for Neg examples\n')
		fout.write('  as One-Class\n')
		fout.write('\n')

		fout.write('Lasso Parameters\n')
		fout.write('method:{}Lasso (standard)\n'.format(textDelim))
		fout.write('alpha:{}{}\n'.format(textDelim, lAlphaPos))
		fout.write('max_iter:{}{}\n'.format(textDelim, lMaxIter))
		fout.write('normalize:{}{}\n'.format(textDelim, lNorm))
		fout.write('positive:{}{}\n'.format(textDelim, 'True'))
		fout.write('fit_intercept:{}{}\n'.format(textDelim, lFitIcpt))
		fout.write('selection:{}{}\n'.format(textDelim, lSelctn))
		fout.write('\n')

		fout.write('Similarity Metric:{}PathSim sum over set\n'.format(textDelim))
		fout.write('Prediction Results\n')
		fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(np.nonzero(cfLasso01.coef_)[0])))
		fout.write('Training score:{}{:3.3f}\n'.format(textDelim, cfLasso01.score(trainSet01, trainLabel01)))
		fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, cfLasso01.score(testSet01, testLabel01)))
		fout.write('\n')
	#end with

#end loop



#TODO: append output file in root of batch ?



print("\nDone.\n")