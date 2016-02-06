# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Predict membership of gene sets using metapaths
#
# This code is to test the idea of using metapath analysis
#	as a one-class classifier. In this script, a portion
#	of a known set will be concealed, and the script will
#	attempt to predict the concealed members.
# Outline:
#	1) 
# ---------------------------------------------------------

import mpLibrary as mp
import time

import numpy as np


# For testing & verification
import random
random.seed(42)



####### ####### ####### ####### 
# PARAMETERS

# Variables/Quantities
percHide = .25		# percent of genes to conceal
nRandSamp = 100		# number of random samples to compare
topKPaths = 10		# number of metapaths to consider
topKGenes = 100		# number of genes to predict


# Input/Output names & locations

useNtwk = 0		# network & samples to use (0 means fake)
if useNtwk == 0 :
#	eName = 'fakeNtwk00_g2e3t10'
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
	sName = 'sample01'
	sPath = 'samplesFake/'
else :
	eName = 'toy2_p3gz'
	ePath = '../networks/'
	sName = 'BERTUCCI_MEDULLARY_VS_DUCTAL_BREAST_CANCER'
	sPath = '../samples/'
#end if

# Path & new directory to save output (& temp files)
oRoot = 'outputFake/'
oRoot = '../output/'
oDirPrefix = 'pred01-' + sName[0:8]

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()
print ""


# 0) Name & create a folder to store output files

oDirectory = mp.nameOutputPath(oRoot, oDirPrefix)
print "Files will be saved to {}".format(oDirectory)
oPath = oRoot + oDirectory



# 1) Load a sample and hide some genes

# Read in the full sample
print "Analyzing metapaths in sample: {}".format(sName)
fullSample = mp.readSampleFiles(sPath + sName, True, True)

#ASIDE: Remove any genes not in the network
inGenes, outGenes = mp.checkGenesInNetwork(ePath,
	eName, fullSample)
print ("Of the {} sample genes,".format(len(fullSample)) +
	" {} are in the network.".format(len(inGenes)) )
del inGenes, outGenes

# Partition into known & concealed sets
gKnown, gHidden = mp.partitionSample(ePath, eName,
	oPath, fullSample, percHide)
print ( "Sample partitioned into {}".format(len(gKnown)) +
	" known and {} concealed genes ...".format(len(gHidden)) )
#print gKnown
#print gHidden
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 2) Identify the paths available

print "Checking what paths are available ..."
pathDict = mp.readKeyFile(ePath, eName)



# 3) Create an array of random samples

# Load the gene-index dict
print "Creating the gene-index dictionary."
geneIndex = mp.readGenesFile(ePath, eName)

# Create N random samples
print ("Choosing {} random samples of".format(nRandSamp) +
	" length {} ...".format(len(gKnown)) )
randSamps = mp.createRandomSamplesArray(nRandSamp,
	len(gKnown), len(geneIndex))
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 4) Calculate the statistics

# Convert the sample into a list of indices
sampIndex = mp.convertToIndices(gKnown, geneIndex)

# Calculate stats for each metapath
print "Calculating statistics ..."
sCount, rMeans, rStDev, zScore, percList = mp.calculateStatistics(
	sampIndex, randSamps, pathDict, ePath, eName)
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 5) Output the metapath stats

# Write to the output file
print "Saving path count & rank data ..."
outFile = mp.writeOutputOneSample(oPath, 'pathlist', eName,
	('known-' + sName), pathDict, sCount, rMeans, rStDev,
	zScore, percList, list([]) )
#print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 6) Choose the top K paths to examine

print "Determining the top {} metapaths ...".format(topKPaths)
bestPaths, bestRanks = mp.chooseTopKPathsSimple(topKPaths, percList, pathDict)

print "Examining the following paths:"
for i in range(len(bestRanks)) :
	print "    {:>3.1f}  {}".format(bestRanks[i], bestPaths[i])
#end loop



# 7) Rank genes by similarity along selected metapaths

# Get the indices for all genes not in the test sample
sampOutdex = [n for n in range(len(geneIndex)) if n not in sampIndex]

# Calculate similarity for each top metapath
simArray = np.empty([len(sampOutdex), len(bestPaths)])
for p in range(len(bestPaths)) :
	simArray[:,p] = mp.applyGroupPathSim(ePath, eName, pathDict[bestPaths[p]], sampIndex)
#end loop

# Write the ranked output to file
outFile = mp.writeItemRanks(oPath, simArray, geneIndex, sampIndex, bestPaths)
print "Saving gene predictions to {}".format(outFile)

	


print "\nDone.\n"