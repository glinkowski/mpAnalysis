# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# A two-step approach to set membership prediction, using
#	a batch approach. This file contains Step One:
#	Rank the paths for multiple samples at once.
#
# The longest part of the prediction process is opening each
#	metapath matrix to compare path counts against random
#	samples. In this approach, all of the samples in the
#	specified folder will be loaded. Then, statistics per
#	path will be calculated in tandem across all samples.
# ---------------------------------------------------------

import mpLibrary as mp
import time
import numpy as np
import gzip


# For testing & verification
#import random
#random.seed(42)



####### ####### ####### ####### 
# PARAMETERS

# Variables/Quantities
percHide = .00		# percent of genes to conceal
nRandSamp = 200		# number of random samples to compare

# options for Node Binning
useBinning = True
nodeBins = [0.3, .66]
binType = 'all'


# Input names & locations
useNtwk = 1		# network & samples to use (0 means fake)
if useNtwk == 0 :
#	eName = 'fakeNtwk00_g2e3t10'
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
	sPath = 'samplesFake/'
else :
	eName = 'all_v1_g2e11t0'
	ePath = '../Dropbox/mp/networks/'
	sPath = '../Dropbox/mp/samples-subset5/'
#	sPath = '../Dropbox/mp/samplesMSIG/'
#end if

# Output path
#oRoot = 'outputFake/'
oRoot = '../Dropbox/mp/output/'
oDirPrefix = 'pred02-batch'


# verbose feedback ?
verbose = False

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()
print ""


#TODO: Does this work? If not ... alternatives?
mp.setParamVerbose(verbose)


# 0) Name & create a folder to store output files
oDirectory = mp.nameOutputPath(oRoot, oDirPrefix)
print "Files will be saved to {}".format(oDirectory)
oPath = oRoot + oDirectory

#TODO: write file w/ experiment parameters


# 1) Load the gene-index dict
print "Creating the gene-index dictionary."
geneDict = mp.readGenesFile(ePath, eName)


# 2) Get list of all samples in folder
sNames = mp.getSampleNamesFromFolder(sPath)
#print sNames


# 3) Read the samples & create random sample sets
oSampLists = list()
rSampArrays = list()

index = 0
for s in sNames :

	# Need a folder for each sample
	oSubDir = oPath+s+'/'

	# Read genes from sample & partition into test/train
	gAll = mp.readSampleFiles(sPath+s, True, True)
	gKnown, gHidden = mp.partitionSample(ePath, eName,
		oSubDir, gAll, percHide)

	print "Analyzing metapaths in sample: {}".format(s)
	print ( "  partitioned into {} known".format(len(gKnown)) +
		" and {} concealed genes ...".format(len(gHidden)) )

	# Convert sample into list of indices
	gIndices = mp.convertToIndices(gKnown, geneDict)
	oSampLists.append(gIndices)

	# Create N random samples
	print ("  choosing {} random samples of".format(nRandSamp) +
		" length {} ...".format(len(gKnown)) )
	if useBinning :
		rSamps = mp.createRandomSamplesBinned(ePath, eName,
			gKnown, geneDict, nodeBins, 'all', nRandSamp)
	else :
		rSamps = mp.createRandomSamplesArray(nRandSamp,
			len(gKnown), len(geneDict))
	#end if
	rSampArrays.append(rSamps)

#TODO: output gene_bins.txt or ... degree_bins ? or ?

#end loop
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)

#print "len of samples list", len(oSampLists)
#print "len of rand samp list", len(rSampArrays)


# 4) Get the list of available paths
print "Checking what paths are available ..."
pathDict = mp.readKeyFile(ePath, eName)
pathList = mp.removeInvertedPaths(pathDict)
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)


# Get expected matrix size
print "Finding the matrix dimensions ..."
fnMatrixZPad = 5
#TODO: pack this into a function
fname = (ePath+eName+"_MetaPaths/" +
	"{}.gz".format(str(0).zfill(fnMatrixZPad)) )
mxSize = 0
with gzip.open(fname, 'rb') as fin :
	for line in fin :
		mxSize += 1
#end with
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)


# 5) Calculate scores 
pathScores = np.zeros([len(pathList), len(sNames)])
for i in xrange(len(pathList)) :

	print "Loading matrix {}, {}, {}".format(i, pathList[i], pathDict[pathList[i]])
	# Load a metapath matrix into memory
	matrix = mp.getPathMatrix(pathDict[pathList[i]], ePath, eName, mxSize)
	print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)

	# Get the percentile/score for this path, per sample
	for j in xrange(len(sNames)) :
		tCount = mp.getPathCountOne(oSampLists[j], matrix)
		tPercent = mp.getPercentile(tCount, rSampArrays[j], matrix)
		pathScores[i,j] = tPercent
	#end loop

	print "  Examined path {}".format(pathList[i])
	print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)
#end loop
#print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)


# 6) output ranked_path.txt per sample
for j in xrange(len(sNames)) :
	oSubDir = oPath+sNames[j]+'/'
	mp.writeRankedPaths(oSubDir, pathScores[:,j], pathDict)
#end loop
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 7) output overall ranked paths for all samples (one file)
textDelim = '\t'
fout = open(oPath+'ranked_paths_all.txt', 'wb')

# Write the file header
fout.write('network:{}{}\n'.format(textDelim, eName))
fout.write('binning:{}'.format(textDelim))
if useBinning :
	fout.write('{}\non edge:{}{}\n'.format(nodeBins, textDelim, binType))
else :
	fout.write('none\n')
#end if
fout.write('\n')

# Write the data header
firstCol = True
for i in xrange(len(sNames)) :
	if not firstCol :
		fout.write(textDelim)
	firstCol = False
	fout.write('{}'.format(sNames[i]))
#end loop
fout.write('\n')

# Write the main data
firstRow = True
for i in xrange(len(pathList)) :
	if not firstRow :
		fout.write('\n')
	firstRow = False

	for j in xrange(len(sNames)) :
		fout.write('{}{}'.format(pathScores[i,j],textDelim))
	#end loop

	fout.write('{}'.format(pathList[i]))
#end loop
fout.close()
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)




print "\nDone.\n"