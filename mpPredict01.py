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



####### ####### ####### ####### 
# PARAMETERS

# The percent of genes to conceal
percHide = .25
# The number of predicted genes to return
numTopK = 100
# How many random samples to examine
numRand = 100

# The network & sample sets to use (0 means fake)
useNtwk = 0

if useNtwk == 0 :
	ename = 'fakeNtwk00_g2e3t10'
	epath = 'networks/'
	sname = 'sample02'
	spath = 'samplesFake/'
else :
	ename = 'toy2_hsa'
	epath = '../networks/'
	sname = 'XXX'
	spath = '../samples/'
#end if

# The directory to save output (& temp)
opath = 'outputFake/pred01-' + sname[0:8] + '/'
#oname = 'pathlist'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()



# 1) Load a sample and hide some genes

# Read in the full sample
print "Analyzing metapaths in sample: {}".format(sname)
fullSample = mp.readSampleFiles(spath + sname, True, True)

#ASIDE: Remove any genes not in the network
inGenes, outGenes = mp.checkGenesInNetwork(epath,
	ename, fullSample)
print ("Of the {} sample genes,".format(len(fullSample)) +
	" {} are in the network.".format(len(inGenes)) )
del inGenes, outGenes

# Partition into known & concealed sets
gKnown, gHidden = mp.partitionSample(epath, ename,
	opath, fullSample, percHide)
#print gKnown
#print gHidden
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 2) Identify the paths available

print "Checking what paths are available ..."
pathDict = mp.readKeyFile(epath, ename)



# 3) Create an array of random samples

# Load the gene-index dict
print "Creating the gene-index dictionary."
geneIndex = mp.readGenesFile(epath, ename)

# Create N random samples
print ("Choosing {} random samples of".format(numRand) +
	" length {} ...".format(len(gKnown)) )
randSamps = mp.createRandomSamplesArray(numRand,
	len(gKnown), len(geneIndex))
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 4) Calculate the statistics

# Convert the sample into a list of indices
sampIndex = mp.convertToIndices(gKnown, geneIndex)

# Calculate stats for each metapath
print "Calculating statistics ..."
sCount, rMeans, rStDev, zScore, percList = mp.calculateStatistics(
	sampIndex, randSamps, pathDict, epath, ename)
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 5) Output the collected data

# Write to the output file
print "Saving path count & rank data ..."
outFile = mp.writeOutputOneSample(opath, 'pathlist', ename,
	('known-' + sname), pathDict, sCount, rMeans, rStDev,
	zScore, percList, list([]) )
#print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)





