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
	eName = 'fakeNtwk00_g2e3t10'
	ePath = 'networks/'
	sName = 'sample02'
	sPath = 'samplesFake/'
else :
	eName = 'toy2_hsa'
	ePath = '../networks/'
	sName = 'XXX'
	sPath = '../samples/'
#end if

# Path & new directory to save output (& temp files)
oRoot = 'outputFake/'
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

## Get the indices for all genes not in the test sample
#sampOutdex = [n for n in range(len(geneIndex))if n not in sampIndex]
##print sampIndex, sampOutdex

#for mp in bestPaths :
#end loop


# test the approach
#path = bestPaths[0]
#colRank = mp.applyGroupPathSim(ePath, eName, pathDict[path], sampIndex)
#print colRank


simArray = np.empty([len(sampIndex), len(bestPaths)])
for p in range(len(bestPaths)) :
	simArray[:,p] = mp.applyGroupPathSim(ePath, eName, pathDict[bestPaths[p]], sampIndex)
#end loop
#print simArray

simSum = np.sum(simArray, axis=0)
simAvg = np.mean(simArray, axis=0)

delim = '\t'
fs = open(oPath + 'scores.txt', 'wb')
fs.write('')
fs.write('Genes scored by similarty to sample ...\n')
fs.write('network:{}{}\n'.format(delim, eName))
fs.write('sample:{}{}\n'.format(delim, sName))
fs.write('known:{}{}\n'.format(delim, len(gKnown)))
fs.write('concealed:{}{}\n'.format(delim, len(gHidden)))
fs.write('\n')

# Get the indices for all genes not in the test sample
sampOutdex = [n for n in range(len(geneIndex)) if n not in sampIndex]
#print sampOutdex
geneList = geneIndex.keys()
geneList.sort()
geneList = [geneList[g] for g in sampOutdex]
del sampOutdex
print geneList

#fs.write('{}'.format(delim))
for bp in bestPaths :
	fs.write('{}{}'.format(delim, bp))
#end loop
fs.write('{}sum{}avg\n'.format(delim, delim))
for i in range( simArray.shape[0] ) :
	fs.write('{}'.format(delim))
	for j in range( simArray.shape[1] ) :
		fs.write( '{}{}'.format(simArray[i,j], delim) )
	#end loop
	fs.write( '{}{}{}{}{}\n'.format(simSum[i], delim,
		simAvg[i], delim, geneList[i]) )
#end loop
fs.write('\n')


	


print "\nDone.\n"