# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Approach 2: learn top paths from PathSim sum, find genes
#	Version 2, Step 1 
#
# The first approach tries to first find the most
#	important metapaths, then rank genes by similarity
#	along those. This approach instead finds the set
#	similarity along each metapath, then tries to learn
#	the most important paths.
#
# Step 1: Build a per-gene feature vector where each entry is
#	the sum of that gene's similarity to each member of
#	the set in question. In this case, the PathSim matrices
#	have already been created as part of the network
#	pre-processing.
# ---------------------------------------------------------

import mpLibrary as mp
import time
import numpy as np
import gzip



####### ####### ####### ####### 
# PARAMETERS

# Variables/Quantities
percHide = [25, 25, 25, 25, 25]  # 5 x 25%
	# percent of genes to conceal


# Input names & locations
useNtwk = 0		# network & samples to use (0 means fake)
if useNtwk == 0 :
#	eName = 'fakeNtwk00_g2e3t10'
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
	sPath = 'samplesFake/'
	oRoot = '../Dropbox/mp/outputFake/'
else :
	eName = 'all_v1_g2e11t0'
	ePath = '../Dropbox/mp/networks/'
	sPath = '../Dropbox/mp/samples-test1/'
	oRoot = '../Dropbox/mp/output/'
#end if

# Output path
oDirPrefix = 'pred03-batch'


# verbose feedback ?
newVerbose = True

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()
print("")

mp.setParamVerbose(newVerbose)



# 0) Name & create a folder to store output files
oDirectory = mp.nameOutputPath(oRoot, oDirPrefix)
print("Files will be saved to {}".format(oDirectory))
oPath = oRoot + oDirectory

# Save experiment parameters to file
fOutput = list()
fOutput.append( ['date', 'network', '% hidden'] )
fOutput.append( [time.strftime('%d/%m/%Y'), eName, percHide] )
fOutputName = 'parameters.gz'
mp.writeGenericLists(oPath, fOutputName, fOutput)



# 1) Load the gene-index dict
print("Creating the gene-index dictionary.")
geneDict = mp.readGenesFile(ePath, eName)



# 2) Get list of all samples in folder
sNames = mp.getSampleNamesFromFolder(sPath)
print sNames






