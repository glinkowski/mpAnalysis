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
lAlpha = 0.05
lMaxIter = 100
lNorm = False
lPos = True


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

#for si in dSubDirs[1:2] :
for si in dSubDirs :

	# Read in the known genes
	print(si)
	gKnown = mp.readFileAsList(si+'known.txt')
#	print(gKnown)
	giKnown = mp.convertToIndices(gKnown, geneDict)
#	print(giKnown)

	# Read in the (modified) group PathSim matrix
	simGroup = mp.readFileAsMatrix(si, fGroupNorm)
	simOrig = mp.readFileAsMatrix(si, fOrigSum)
#	print(simGroup.shape, simOrig.shape)
#	print(np.amax(simGroup,axis=1)[0:12])
#	print(np.amax(simOrig,axis=1)[0:12])

	# Extract the positive data on which to perform LASSO
	posGroup = simGroup[giKnown,:]
	posOrig = simOrig[giKnown,:]
	posLabel = np.ones(len(giKnown))
#	print(posGroup[0:3,0:10])

	# Select & extract random counter-examples
	nExamples = min( len(giKnown), (len(geneDict) - len(giKnown)) )
	giUnknown = [g for g in geneDict.values() if g not in giKnown]
	giNegTrain = random.sample(giUnknown, nExamples)
	negGroup = simGroup[giUnknown,:]
	negOrig = simOrig[giUnknown,:]
	negLabel = np.zeros(len(giUnknown))
#	print(negGroup[0:3,0:10])

	print("  ... performing regression ...")
#	print posLabel.shape, negLabel.shape
	# Combine into the full training set
	trainGroup = np.vstack( (posGroup, negGroup) )
#	print trainGroup.shape
	trainOrig = np.vstack( (posOrig, negOrig) )

	trainLabel = np.vstack( (np.reshape(posLabel, (len(posLabel), 1)),
		np.reshape(negLabel, (len(negLabel), 1))) )

	# Perform LASSO (x2)
	lGroup = lm.Lasso(alpha=lAlpha, max_iter=lMaxIter,
		normalize=lNorm, positive=lPos)
	lGroup.fit(trainGroup, trainLabel)
#	print(lGroup.coef_.shape, lGroup.coef_)

	lOrig = lm.Lasso(alpha=lAlpha, max_iter=lMaxIter,
		normalize=lNorm, positive=lPos)
	lOrig.fit(trainOrig, trainLabel)
#	print(lOrig.coef_.shape, lOrig.coef_)

	# Extract the weights and corresponding paths
#	nodeDT = np.dtype('a30')

#	print lGroup.coef_

	iCoefGroup = np.nonzero(lGroup.coef_)[0]
#	print(iCoefGroup[0])
#	pGroup = np.recarray( len(iCoefGroup), dtype=[('path', nodeDT), ('weight', 'f4')] )
	pGroup = np.recarray( len(iCoefGroup), dtype=[('path', 'i4'), ('weight', 'f4')] )
	row = 0
	for c in iCoefGroup :
#		print(c)
#		print(pathNames[c], lGroup.coef_[c])
#		pGroup[row] = (pathNames[c], lGroup.coef_[c])
		pGroup[row] = (c, lGroup.coef_[c])
		row += 1
	pGroup[::-1].sort(order=['weight', 'path'])	# sort by descending wieght

	iCoefOrig = np.nonzero(lOrig.coef_)[0]
	pOrig = np.recarray( len(iCoefOrig), dtype=[('path', 'i4'), ('weight', 'f4')] )
	row = 0
	for c in iCoefOrig :
#		pOrig[row] = (pathNames[c], lOrig.coef_[c])
		pOrig[row] = (c, lOrig.coef_[c])
		row += 1
	pOrig[::-1].sort(order=['weight', 'path'])	# sort by descending wieght


	# Output the coefficients (x2)
	textDelim = '\t'

	fPrefix = 'top_paths_Lasso-'+fGroupNorm.rstrip('.txtgz')
	fname = mp.nameOutputFile(si, fPrefix)
	print("Saving top paths to file {}".format(fname))
#	print("  in directory {}".format(si))
	with open(si+fname, 'wb') as fout :
	#	fout = open(si+fname, 'wb')
		fout.write('alpha:{0}{1}{0}max_iter:{0}{2}{0}'.format(textDelim, lAlpha, lMaxIter) +
			'normalize:{0}{1}{0}positive:{0}{2}{0}'.format(textDelim, lNorm, lPos))
		for row in xrange(len(pGroup)) :
			fout.write('\n{}{}{}'.format(pGroup['weight'][row],
				textDelim, pathNames[pGroup['path'][row]]))
	#		fout.write('\n{}{}{}'.format(row[1], textDelim, row[0]))
	#		fout.write('\n{}{}{}'.format(row['path'], textDelim, row['weight']))
	#		fouta.write("{}{}{}".format(rankList['score'][i], textDelim, rankList['names'][i]))
	#end with
#	fout.close()

	fPrefix = 'top_paths_Lasso-'+fOrigSum.rstrip('.txtgz')
	fname = mp.nameOutputFile(si, fPrefix)
	print("Saving top paths to file {}".format(fname))
#	print("  in directory {}".format(si))
	with open(si+fname, 'wb') as fout :
		fout.write('alpha:{0}{1}{0}max_iter:{0}{2}{0}'.format(textDelim, lAlpha, lMaxIter) +
			'normalize:{0}{1}{0}positive:{0}{2}{0}'.format(textDelim, lNorm, lPos))
		for row in xrange(len(pOrig)) :
			fout.write('\n{}{}{}'.format(pOrig['weight'][row],
				textDelim, pathNames[pOrig['path'][row]]))
	#end with
#end loop


# rank genes

# compare the two metrics

