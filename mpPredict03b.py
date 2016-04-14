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
	dDir = 'pred03-batch-001'
#end if


# File names for similarity metrics
#fGroupNorm = 'Pxy.gz'
fGroupNorm = 'Pxy-mod-norm.gz'
fOrigSum = 'SxySum.gz'


# LASSO params
lAlpha = 0.09
lMaxIter = 10000
lNorm = True
lPos = True
lFitIcpt = True
lSelctn = 'random' # random vs cyclic


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

for si in dSubDirs[4:6] :
#for si in dSubDirs :

	print("\n", si)


	# ####### ####### ####### #######
	# Prepare the test/train vectors & labels
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
	#	(as one-class) using rand sample of Unknown
	nExamples = min( 4 * len(giKnown), (len(geneDict) - len(giKnown)) )
	giRandNeg = random.sample(giUnknown, nExamples)
#	giRandNeg = random.sample(giTrueNeg, nExamples)

	negTrainG = vectG[giRandNeg]
	negTestG = vectG[giTrueNeg]

	negTrainO = vectO[giRandNeg]
	negTestO = vectO[giTrueNeg]

#	negTrainLabel = np.zeros(len(giRandNeg))
	negTrainLabel = np.ones( (len(giRandNeg), 1) ) * -100
#	negTrainLabel = np.ones(len(giRandNeg)) * -1

#	negTestLabel = np.zeros(len(giTrueNeg))
	negTestLabel = np.ones( (len(giTrueNeg), 1) ) * -100
#	negTestLabel = np.ones(len(giTrueNeg)) * -1

	# Combine to create the full train & test data sets
	trainG = np.vstack( (posTrainG, negTrainG) )
	trainO = np.vstack( (posTrainO, negTrainO) )
	trainLabel = np.vstack( (posTrainLabel, negTrainLabel) )

	testG = np.vstack( (posTestG, negTestG) )
	testO = np.vstack( (posTestO, negTestO) )
	testLabel = np.vstack( (posTestLabel, negTestLabel) )


	# ####### ####### ####### #######
	# Perform the regression analysis
	print("  ... performing regression ...")

	# Train LASSO
#TODO: exp w/ diff types: LassoCV, ElasticNet, MultiTaskElasticNet, MultiTaskLasso ..?
	lassoG = lm.Lasso(alpha=lAlpha, max_iter=lMaxIter, normalize=lNorm,
		positive=lPos, fit_intercept=lFitIcpt, selection=lSelctn)
	lassoG.fit(trainG, trainLabel)


	# How well does it work ??
	# Get score from training data
	scoreG = lassoG.score(trainG, trainLabel)
	print("On training data {}, score: {}".format(fGroupNorm, scoreG))
	print("  using {} coefficients".format( len(np.nonzero(lassoG.coef_)[0]) ))


	# On the full set ?? ...
	scoreG = lassoG.score(testG, testLabel)
	print("On test data {}, score: {}".format(fGroupNorm, scoreG))
	predLabel = lassoG.predict(testG)
	print("Some labels for the hidden set:")
#	print(predLabel[0:len(posTestLabel)])
	print(predLabel[0:5])
	print("Some labels for the true neg set:")
	print(predLabel[len(posTestLabel):(len(posTestLabel) + 5)])





	iCoefGroup = np.nonzero(lassoG.coef_)[0]
#	print(iCoefGroup[0])
#	pGroup = np.recarray( len(iCoefGroup), dtype=[('path', nodeDT), ('weight', 'f4')] )
	pGroup = np.recarray( len(iCoefGroup), dtype=[('path', 'i4'), ('weight', 'f4')] )
	row = 0
	for c in iCoefGroup :
#		print(c)
#		print(pathNames[c], lGroup.coef_[c])
#		pGroup[row] = (pathNames[c], lGroup.coef_[c])
		pGroup[row] = (c, lassoG.coef_[c])
		row += 1
	pGroup[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

#	iCoefOrig = np.nonzero(lOrig.coef_)[0]
#	pOrig = np.recarray( len(iCoefOrig), dtype=[('path', 'i4'), ('weight', 'f4')] )
#	row = 0
#	for c in iCoefOrig :
##		pOrig[row] = (pathNames[c], lOrig.coef_[c])
#		pOrig[row] = (c, lOrig.coef_[c])
#		row += 1
#	pOrig[::-1].sort(order=['weight', 'path'])	# sort by descending wieght


	# Output the coefficients (x2)
	textDelim = '\t'

	fPrefix = 'top_paths_Lasso-'+fGroupNorm.rstrip('.txtgz')
	fname = mp.nameOutputFile(si, fPrefix)
	print("Saving top paths to file {}".format(fname))
#	print("  in directory {}".format(si))
	with open(si+fname, 'wb') as fout :
	#	fout = open(si+fname, 'wb')
		fout.write('alpha:{0}{1}{0}max_iter:{0}{2}{0}'.format(textDelim, lAlpha, lMaxIter) +
			'normalize:{0}{1}{0}positive:{0}{2}{0}'.format(textDelim, lNorm, lPos) +
			'fit_intercept:{0}{1}{0}selection:{0}{2}{0}'.format(textDelim, lFitIcpt, lSelctn))
		for row in xrange(len(pGroup)) :
			fout.write('\n{}{}{}'.format(pGroup['weight'][row],
				textDelim, pathNames[pGroup['path'][row]]))
	#		fout.write('\n{}{}{}'.format(row[1], textDelim, row[0]))
	#		fout.write('\n{}{}{}'.format(row['path'], textDelim, row['weight']))
	#		fouta.write("{}{}{}".format(rankList['score'][i], textDelim, rankList['names'][i]))
	#end with
#	fout.close()

#	fPrefix = 'top_paths_Lasso-'+fOrigSum.rstrip('.txtgz')
#	fname = mp.nameOutputFile(si, fPrefix)
#	print("Saving top paths to file {}".format(fname))
##	print("  in directory {}".format(si))
#	with open(si+fname, 'wb') as fout :
#		fout.write('alpha:{0}{1}{0}max_iter:{0}{2}{0}'.format(textDelim, lAlpha, lMaxIter) +
#			'normalize:{0}{1}{0}positive:{0}{2}{0}'.format(textDelim, lNorm, lPos))
#		for row in xrange(len(pOrig)) :
#			fout.write('\n{}{}{}'.format(pOrig['weight'][row],
#				textDelim, pathNames[pOrig['path'][row]]))
#	#end with




#end loop


# rank genes

# compare the two metrics

