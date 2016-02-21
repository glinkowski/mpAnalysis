# ---------------------------------------------------------
# author: Greg Linkowski
# node binning by Shengliang Dai
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# Find all the metapaths present in a specific sample
#
# This is a basic demo. Load a sample file and rank all the
#   paths found between the genes in that sample. Do this
#   by comparing against random samples of the same size to
#   see which paths are more pronounced in the sample.
# Outline:
#	1) Read in the sample file
#	2) Identify the path matricies from the network
#	3) Create a number of random samples
#	4) Calculate the statistics for this sample
#	5) Output the results to a file
# ---------------------------------------------------------

import mpFindFuncs as ff
import mpLibrary as mpl
#import mpFindFuncs02 as ff2
import time

## REMOVE: The following is for verification.
#import random
#random.seed(42)



####### ####### ####### ####### 
# PARAMETERS

# The network to use and directory path
# The sample to test and path
useRealData = True

if not useRealData :
	ename = 'fakeNtwk00_g2e3t10'
	ename = 'fakeNtwk01_g3e4t1'
	epath = 'networks/'
	sname = 'sample04'
	spath = 'samplesFake/'
else :
	ename = 'all_v1_g2e11t0'
#	ename = 'toy2_p3gz'
	epath = '../networks/'
	sname = 'CAMPS_COLON_CANCER_COPY_NUMBER'
	spath = '../samplesMSIG/'
#end if


# Where to store the output
oname = 'mpf02-' + ename + '-' + sname
opath = '../Dropbox/mp/output/'


# How many random samples to examine
numRand = 100

# NODE BINNING --
# Threshold percentages for binning the genes
binThresholds = [.33, .66]

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION
print ""

tstart = time.time()



# 1) Read in the sample file

print "Finding metapaths in sample: {}".format(sname)
sampGenes = ff.readSampleFiles(spath + sname, True, True)
#print sampGenes
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 2) Identify the paths available

print "Checking what paths are available ..."
pathDict = ff.readKeyFile(epath, ename)
#print mpDict
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 3) Create an array of random samples

# Check which genes are actually in the network
inGenes, outGenes = ff.checkGenesInNetwork(epath, ename, sampGenes)
print ("Of the {} sample genes,".format(len(sampGenes)) +
	" {} are in the network.".format(len(inGenes)) )
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)

# Load the gene-index dict
print "Creating the gene-index dictionary."
geneIndex = ff.readGenesFile(epath, ename)
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)

# Create N random samples
print ("Choosing {} random samples of".format(numRand) +
	" length {} ...".format(len(inGenes)) )

######## ######## ######## ######## 
#TODO: Replace the following line.
#	Instead, use node binning to select the random samples.
#randSamps = ff2.sampleMatrix(epath, ename, sampGenes, spath, sname)
######## ######## ######## ######## 

#randSamps = ff.createRandomSamplesArray(numRand,
#	len(inGenes), len(geneIndex))
randSamps = mpl.createRandomSamplesBinned(epath, ename, inGenes,
	geneIndex, binThresholds, 'all', numRand)



print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 4) Calculate the statistics

# Convert the sample into a list of indices
sampIndex = ff.convertToIndices(inGenes, geneIndex)
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)

# Calculate stats for each metapath
print "Calculating statistics ..."
sCount, rMeans, rStDev, zScore, percList = ff.calculateStatistics(
	sampIndex, randSamps, pathDict, epath, ename)
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 5) Output the collected data

# Write to the output file
print "Saving data to ..."
outFile = ff.writeOutputOneSample(opath, oname, ename, sname,
	pathDict, sCount, rMeans, rStDev, zScore, percList, outGenes)
print "    ... {}".format(outFile)
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



print "\nDone.\n"
