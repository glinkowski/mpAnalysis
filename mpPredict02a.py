# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Approach 1: path counting & percentile beat
#	Version 2, Step 1
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
percHide = .25		# percent of genes to conceal
nRandSamp = 100		# number of random samples to compare

# options for Node Binning
useBinning = True
nodeBins = [0.4, .85]
binType = 'all'


# Input names & locations
useNtwk = 1		# network & samples to use (0 means fake)
if useNtwk == 0 :
#	eName = 'fakeNtwk00_g2e3t10'
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
	sPath = 'samplesFake/'
	oRoot = '../Dropbox/mp/outputFake/'
else :
	eName = 'all_v3beta_g2e9t0'
	ePath = '../Dropbox/mp/networks/'
	sPath = '../Dropbox/mp/samples-4subs/subset02/'
#	sPath = '../Dropbox/mp/samplesMSIG/'
	oRoot = '../Dropbox/mp/output/'
#end if

# Output path
oDirPrefix = 'pred02a-batch'


# verbose feedback ?
verbose = False

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()
print("")

mp.setParamVerbose(verbose)



# 0) Name & create a folder to store output files
oDirectory = mp.nameOutputPath(oRoot, oDirPrefix)
print("Files will be saved to {}".format(oDirectory))
oPath = oRoot + oDirectory



# 1) Load the gene-index dict
print("Creating the gene-index dictionary.")
geneDict = mp.readGenesFile(ePath, eName)



# 2) Get list of all samples in folder
sNames = mp.getSampleNamesFromFolder(sPath)
#printsNames



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

	print("Analyzing metapaths in sample: {}".format(s))
	print( "  partitioned into {} known".format(len(gKnown)) +
		" and {} concealed genes ...".format(len(gHidden)) )

	# Convert sample into list of indices
	gIndices = mp.convertToIndices(gKnown, geneDict)
	oSampLists.append(gIndices)

	# Create N random samples
	print("  choosing {} random samples of".format(nRandSamp) +
		" length {} ...".format(len(gKnown)) )
	if useBinning :
		rSamps = mp.createRandomSamplesBinned(ePath, eName,
			gKnown, geneDict, nodeBins, 'all', nRandSamp)
	else :
		rSamps = mp.createRandomSamplesArray(nRandSamp,
			len(gKnown), len(geneDict))
	#end if
	rSampArrays.append(rSamps)
#end loop
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))



# 4) Get the list of available paths
print("Checking what paths are available ...")
pathDict = mp.readKeyFile(ePath, eName)
pathList = mp.removeInvertedPaths(pathDict)
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))


# Get expected matrix size
print("Finding the matrix dimensions ...")
fnMatrixZPad = 5
#TODO: pack this into a function
fname = (ePath+eName+"_MetaPaths/" +
	"{}.gz".format(str(0).zfill(fnMatrixZPad)) )
mxSize = 0
with gzip.open(fname, 'rb') as fin :
	for line in fin :
		mxSize += 1
#end with
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))



# 5) Calculate scores 
pathScores = np.zeros([len(pathList), len(sNames)])
pathScores2 = np.zeros([len(pathList), len(sNames)])
for i in range(len(pathList)) :

	print("Loading matrix {}, {}, {}".format(i, pathList[i], pathDict[pathList[i]]))
	# Load a metapath matrix into memory
	matrix = mp.getPathMatrix(pathDict[pathList[i]], ePath, eName, mxSize)
	print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))

	# Get the percentile/score for this path, per sample
	for j in range(len(sNames)) :
		tCount = mp.getPathCountOne(oSampLists[j], matrix)
		tPercent, tPercDiff = mp.getPercentile(tCount, rSampArrays[j], matrix)
		pathScores[i,j] = tPercent
#		tPercDiff = mp.getPercentDifference(tCount, rSampArrays[j], matrix)
		pathScores2[i,j] = tPercDiff
	#end loop

	print("  Examined path {}".format(pathList[i]))
	print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))
#end loop
#print"    --elapsed time: {:.3} (s)".format(time.time()-tstart)


#TODO: save pathScores2

# 6) output ranked_path.txt per sample
for j in range(len(sNames)) :
	oSubDir = oPath+sNames[j]+'/'
	mp.writeRankedPaths(oSubDir, 'percent', pathScores[:,j], pathDict)
	mp.writeRankedPaths(oSubDir, 'difference', pathScores2[:,j], pathDict)
#end loop
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))



# 7) output overall ranked paths for all samples (one file)
textDelim = '\t'
fout = open(oPath+'ranked_paths_all-percent.txt', 'w')

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
for i in range(len(sNames)) :
	if not firstCol :
		fout.write(textDelim)
	firstCol = False
	fout.write('{}'.format(sNames[i]))
#end loop
fout.write('\n')

# Write the main data
firstRow = True
for i in range(len(pathList)) :
	if not firstRow :
		fout.write('\n')
	firstRow = False

	for j in range(len(sNames)) :
		fout.write('{}{}'.format(pathScores[i,j],textDelim))
	#end loop

	fout.write('{}'.format(pathList[i]))
#end loop
fout.close()
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))


fout = open(oPath+'ranked_paths_all-difference.txt', 'w')
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
for i in range(len(sNames)) :
	if not firstCol :
		fout.write(textDelim)
	firstCol = False
	fout.write('{}'.format(sNames[i]))
#end loop
fout.write('\n')
# Write the main data
firstRow = True
for i in range(len(pathList)) :
	if not firstRow :
		fout.write('\n')
	firstRow = False
	for j in range(len(sNames)) :
		fout.write('{}{}'.format(pathScores2[i,j],textDelim))
	#end loop
	fout.write('{}'.format(pathList[i]))
#end loop
fout.close()


# 8) write parameters, network, sample names to file
with open(oPath + 'parameters.txt', 'w') as fout :
	fout.write('network\t{}\n'.format(eName))
	fout.write('ntwk path\t{}\n'.format(ePath))
	fout.write('\n')
	fout.write('percent hidden\t{:2.1%}\n'.format(percHide))
	fout.write('random samples\t{}\n'.format(nRandSamp))
	fout.write('used node binning\t{}\n'.format(useBinning))
	fout.write('bin cutoffs\t{}\n'.format(nodeBins))
	fout.write('bin edge type\t{}\n'.format(binType))
	fout.write('\n')
	fout.write('runtime (s)\t{:1.3}\n'.format(time.time()-tstart))
	fout.write('samples (#)\t{}\n'.format(len(sNames)))
	fout.write('Samples in this batch:')
	for sn in sNames :
		fout.write('\n\t{}'.format(sn))
#end with



print("\nDone.\n")