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
	dDir = 'pred03-batch-002'
#end if

# whether to write output files
writeOutput = True


# File names for similarity metrics
#fGroupNorm = 'Pxy.gz'
fGroupNorm = 'Pxy-mod-norm.gz'
fOrigSum = 'SxySum.gz'


# how to sample negative train/test set
sampleAsOneClass = True


# LASSO params
if sampleAsOneClass :
	lAlpha = 0.05	# good for asOneClass
else :
	lAlpha = 0.005 	# find good TwoClass alpha
#end if
lMaxIter = 10000
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

	print("\n{}".format(si))



	# ####### ####### ####### #######
	# 3a) Prepare the test/train vectors & labels
	print("  ... creating train & test data ...")

	# Read in the (modified) group PathSim matrix
	vectG = mp.readFileAsMatrix(si, fGroupNorm)
	vectO = mp.readFileAsMatrix(si, fOrigSum)

	# Create index lists for Known, Hidden, Unknown, TrueNeg
	gKnown = mp.readFileAsList(si+'known.txt')
	giKnown = mp.convertToIndices(gKnown, geneDict)
	gHidden = mp.readFileAsList(si+'concealed.txt')
	giHidden = mp.convertToIndices(gHidden, geneDict)
	giUnknown = [g for g in geneDict.values() if g not in giKnown]
	giTrueNeg = [g for g in giUnknown if g not in giHidden]
	print len(gKnown), len(gHidden), len(giUnknown), len(giTrueNeg)

	# Extract the vectors for the pos sets
	posTrainG = vectG[giKnown,:]
	posTestG = vectG[giHidden,:]

	posTrainO = vectO[giKnown,:]
	posTestO = vectO[giHidden,:]

	posTrainLabel = np.ones( (len(giKnown), 1) ) * 100
#	posTrainLabel = np.ones(len(giKnown))

	posTestLabel = np.ones( (len(giHidden), 1) ) * 100
#	posTestLabel = np.ones(len(giHidden))


	# Extract the vectors for the neg sets
#TODO: best way to handle/produce negative train data ??
	if sampleAsOneClass == True :
		# as one-class: train with rand samp from Unknown
		#		test with all Unknown (TrueNeg + Hidden)
		nExamples = min( 2 * len(giKnown), (len(geneDict) - len(giKnown)) )
		giTrainNeg = random.sample(giUnknown, nExamples)

		negTrainG = vectG[giTrainNeg]
		negTrainO = vectO[giTrainNeg]
		negTrainLabel = np.ones( (len(giTrainNeg), 1) ) * -100

		negTestG = vectG[giTrueNeg]
		negTestO = vectO[giTrueNeg]
		negTestLabel = np.ones( (len(giTrueNeg), 1) ) * -100

		giTest = np.hstack( (giHidden, giTrueNeg) )
	else :
		# as two-class: train rand samp from TrueNeg (1/2)
		#		test with remaining (TrueNeg * 1/2 + Hidden)
		nExamples = int( len(giTrueNeg) / 2 )
		giTrainNeg = random.sample(giTrueNeg, nExamples)
		giTestNeg = [g for g in giTrueNeg if g not in giTrainNeg]

		negTrainG = vectG[giTrainNeg]
		negTrainO = vectO[giTrainNeg]
		negTrainLabel = np.ones( (len(giTrainNeg), 1) ) * -100

		negTestG = vectG[giTestNeg]
		negTestO = vectO[giTestNeg]
		negTestLabel = np.ones( (len(giTestNeg), 1) ) * -100

		giTest = np.hstack( (giHidden, giTestNeg) )
	#end if


	# Combine to create the full train & test data sets
	trainG = np.vstack( (posTrainG, negTrainG) )
	trainO = np.vstack( (posTrainO, negTrainO) )
	trainLabel = np.vstack( (posTrainLabel, negTrainLabel) )

	testG = np.vstack( (posTestG, negTestG) )
	testO = np.vstack( (posTestO, negTestO) )
	testLabel = np.vstack( (posTestLabel, negTestLabel) )



	# ####### ####### ####### #######
	# 3b) Perform the regression analysis
	print("  ... performing regression ...")

	# Train LASSO	

	# Some versions want the labels reshaped
	trainLabel = np.reshape(trainLabel, [trainLabel.shape[0],])

#TODO: exp w/ diff types: LassoCV, ElasticNet, MultiTaskElasticNet, MultiTaskLasso ..?
	lassoG = lm.Lasso(alpha=lAlpha, max_iter=lMaxIter, normalize=lNorm,
	 	positive=lPos, fit_intercept=lFitIcpt, selection=lSelctn)
#	lassoG = lm.LassoCV(max_iter=lMaxIter, normalize=lNorm,
#		positive=lPos, fit_intercept=lFitIcpt, selection=lSelctn)
#	lassoG = lm.ElasticNetCV(max_iter=lMaxIter, normalize=lNorm,
#	 	positive=lPos, fit_intercept=lFitIcpt)
	lassoG.fit(trainG, trainLabel)


	# How well does it work ??
	# Get score from training data
	scoreG = lassoG.score(trainG, trainLabel)
	print("On training data {}, score: {}".format(fGroupNorm, scoreG))
	print("  using {} coefficients".format( len(np.nonzero(lassoG.coef_)[0]) ))


	# On the full set ?? ...
	scoreG = lassoG.score(testG, testLabel)
	print("On test data {}, score: {}".format(fGroupNorm, scoreG))
	predLabelG = lassoG.predict(testG)
#	print("Some labels for the hidden set:")
##	print(predLabel[0:len(posTestLabel)])
#	print(predLabel[0:5])
#	print("Some labels for the true neg set:")
#	print(predLabel[len(posTestLabel):(len(posTestLabel) + 5)])


	# Train LASSO (2nd metric)
#TODO: exp w/ diff types: LassoCV, ElasticNet, MultiTaskElasticNet, MultiTaskLasso ..?
	lassoO = lm.Lasso(alpha=lAlpha, max_iter=lMaxIter, normalize=lNorm,
		positive=lPos, fit_intercept=lFitIcpt, selection=lSelctn)
#	lassoO = lm.LassoCV(max_iter=lMaxIter, normalize=lNorm,
#		positive=lPos, fit_intercept=lFitIcpt, selection=lSelctn)
#	lassoO = lm.ElasticNetCV(max_iter=lMaxIter, normalize=lNorm,
#		positive=lPos, fit_intercept=lFitIcpt)
	lassoO.fit(trainO, trainLabel)


	# How well does it work ??
	# Oet score from training data
	scoreO = lassoO.score(trainO, trainLabel)
	print("On training data {}, score: {}".format(fOrigSum, scoreO))
	print("  using {} coefficients".format( len(np.nonzero(lassoO.coef_)[0]) ))


	# On the full set ?? ...
	scoreO = lassoO.score(testO, testLabel)
	print("On test data {}, score: {}".format(fOrigSum, scoreO))
	predLabelO = lassoO.predict(testO)
#	print("Some labels for the hidden set:")
##	print(predLabel[0:len(posTestLabel)])
#	print(predLabel[0:5])
#	print("Some labels for the true neg set:")
#	print(predLabel[len(posTestLabel):(len(posTestLabel) + 5)])



	# ####### ####### ####### #######
	# 3c) Output selected path results to file

	# Save the selected paths & scores/weights
	if writeOutput :

		textDelim = '\t'


		# Output the feature coefficients (mp weights)
		#	for the first metric
		iCoefGroup = np.nonzero(lassoG.coef_)[0]
	#	print(iCoefGroup[0])
		pGroup = np.recarray( len(iCoefGroup), dtype=[('path', 'i4'), ('weight', 'f4')] )
		row = 0
		for c in iCoefGroup :
			pGroup[row] = (c, lassoG.coef_[c])
			row += 1
		pGroup[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

		# write the file
		fPrefix = 'top_paths_Lasso-'+fGroupNorm.rstrip('.txtgz')
		fname = mp.nameOutputFile(si, fPrefix)
		print("Saving top paths to file {}".format(fname))
	#	print("  in directory {}".format(si))
		with open(si+fname, 'wb') as fout :
			fout.write('alpha:{0}{1}{0}max_iter:{0}{2}{0}'.format(textDelim, lAlpha, lMaxIter) +
				'normalize:{0}{1}{0}positive:{0}{2}{0}'.format(textDelim, lNorm, lPos) +
				'fit_intercept:{0}{1}{0}selection:{0}{2}{0}'.format(textDelim, lFitIcpt, lSelctn))
			for row in xrange(len(pGroup)) :
				fout.write('\n{}{}{}'.format(pGroup['weight'][row],
					textDelim, pathNames[pGroup['path'][row]]))
		#end with

		# Output the feature coefficients (mp weights)
		#	for the second metric
		iCoefOrig = np.nonzero(lassoO.coef_)[0]
		pOrig = np.recarray( len(iCoefOrig), dtype=[('path', 'i4'), ('weight', 'f4')] )
		row = 0
		for c in iCoefOrig :
			pOrig[row] = (c, lassoO.coef_[c])
			row += 1
		pOrig[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

		# write the file
		fPrefix = 'top_paths_Lasso-'+fOrigSum.rstrip('.txtgz')
		fname = mp.nameOutputFile(si, fPrefix)
		print("Saving top paths to file {}".format(fname))
	#	print("  in directory {}".format(si))
		with open(si+fname, 'wb') as fout :
			fout.write('alpha:{0}{1}{0}max_iter:{0}{2}{0}'.format(textDelim, lAlpha, lMaxIter) +
				'normalize:{0}{1}{0}positive:{0}{2}{0}'.format(textDelim, lNorm, lPos) +
				'fit_intercept:{0}{1}{0}selection:{0}{2}{0}'.format(textDelim, lFitIcpt, lSelctn))
			for row in xrange(len(pOrig)) :
				fout.write('\n{}{}{}'.format(pOrig['weight'][row],
					textDelim, pathNames[pOrig['path'][row]]))
		#end with
	#end if



	# ####### ####### ####### #######
	# 3d) Output ranked gene results to file

	# Save the selected paths & scores/weights
	if writeOutput :

		textDelim = '\t'

		geneList = geneDict.keys()
		geneList.sort()

		# Sort the genes by (inverse) rank
		#	first metric
		rGenesG = np.recarray( len(predLabelG), dtype=[('gene', 'i4'), ('rank', 'f4')] )
		rGenesG['gene'] = giTest
		rGenesG['rank'] = predLabelG
		rGenesG[::-1].sort(order=['rank','gene'])

		# write the file
		fPrefix = 'ranked_genes-'+fGroupNorm.rstrip('.txtgz')
		fname = mp.nameOutputFile(si, fPrefix)
		print("Saving ranked genes to file {}".format(fname))
		with open(si+fname, 'wb') as fout :
			firstRow = True
			for row in xrange(len(rGenesG)) :
				if not firstRow :
					fout.write('\n')
				fout.write('{:3.3f}{}{}'.format(rGenesG['rank'][row],
					textDelim, geneList[rGenesG['gene'][row]]))
				firstRow = False
		#end with

		# Sort the genes by (inverse) rank
		#	second metric
		rGenesO = np.recarray( len(predLabelO), dtype=[('gene', 'i4'), ('rank', 'f4')] )
		rGenesO['gene'] = giTest
		rGenesO['rank'] = predLabelO
		rGenesO[::-1].sort(order=['rank','gene'])

		# write the file
		fPrefix = 'ranked_genes-'+fOrigSum.rstrip('.txtgz')
		fname = mp.nameOutputFile(si, fPrefix)
		print("Saving ranked genes to file {}".format(fname))
		with open(si+fname, 'wb') as fout :
			firstRow = True
			for row in xrange(len(rGenesO)) :
				if not firstRow :
					fout.write('\n')
				fout.write('{:3.3f}{}{}'.format(rGenesO['rank'][row],
					textDelim, geneList[rGenesO['gene'][row]]))
				firstRow = False
		#end with
	#end if



	# ####### ####### ####### #######
	# 3e) Output parameters & scores to file

	# Save the parameters & results
	if writeOutput :

		textDelim = '\t'

#		fPrefix = 'parameters-'+fGroupNorm.rstrip('.txtgz')
#		fname = mp.nameOutputFile(si, fPrefix)
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
			fout.write('alpha:{}{}\n'.format(textDelim, lAlpha))
			fout.write('max_iter:{}{}\n'.format(textDelim, lMaxIter))
			fout.write('normalize:{}{}\n'.format(textDelim, lNorm))
			fout.write('positive:{}{}\n'.format(textDelim, lPos))
			fout.write('fit_intercept:{}{}\n'.format(textDelim, lFitIcpt))
			fout.write('selection:{}{}\n'.format(textDelim, lSelctn))
			fout.write('\n')

			fout.write('Similarity Metric:{}{}\n'.format(textDelim, fGroupNorm))
			fout.write('Prediction Results\n')
			fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(np.nonzero(lassoG.coef_)[0])))
			fout.write('Training score:{}{:3.3f}\n'.format(textDelim, lassoG.score(trainG, trainLabel)))
			fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, lassoG.score(testG, testLabel)))
			fout.write('\n')

			fout.write('Similarity Metric:{}{}\n'.format(textDelim, fGroupNorm))
			fout.write('Prediction Results\n')
			fout.write('nonzero coefficients:{}{}\n'.format(textDelim, len(np.nonzero(lassoO.coef_)[0])))
			fout.write('Training score:{}{:3.3f}\n'.format(textDelim, lassoO.score(trainO, trainLabel)))
			fout.write('Testing score:{}{:3.3f}\n'.format(textDelim, lassoO.score(testO, testLabel)))
			fout.write('\n')
		#end with
	#end if

#end loop



print("\nDone.\n")