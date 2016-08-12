# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# Pre-Processing of the network
#	Create features from the non-gene network nodes (terms)
#
# This assumes the pre-processing has already been run.
#	From the normalized edge file, create a gene-by-term
#	matrix where each entry is the edge weight indicating
#	how strongly that gene is defined as belonging to
#	that term.
# ---------------------------------------------------------

import preProcFuncs as pp
import matplotlib.pyplot as plt
import numpy as np
import gzip
import time
import mpLibrary as mp




####### ####### ####### ####### 
# PARAMETERS

# The network to use and directory path
useFakeData = False
if useFakeData :
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
else :
	eName = 'all_v3beta_g2e9t0'
	ePath = '../Dropbox/mp/networks/'
#end if

# Wheather to save uncompressed text copy of matrix
saveAsText = False

# The data type to use for the feature matrix
featDType = np.float32

####### ####### ####### ####### 




####### ####### ####### ####### 
# ANCILLARY FUNCTIONS

def reverseName(pName) :
	pv = pName.split('-')
	rv = pv[::-1]
	revName = ''
	for part in rv :
		revName = revName + part + '-'
#	revName.rstrip('-')
	revName = revName[0:-1]

	return revName
#end def ############### 
def concatenatePaths(path01, path02) :
	if not path01.endswith('/') :
		path01 = path01 + '/'
	if not path02.endswith('/') :
		path02 = path02 + '/'

	return (path01 + path02)
#end def ############### 




####### ####### ####### ####### 
# PRIMARY FUNCTION
# Function: create feature matrix of non-gene terms
# Input ----
#	eList, array: normalized network edge list
#	indirect, str list: names of edge types with non-gene nodes
#	gList, str list: names of the genes in the network
# Returns ----
def createFeatureTerms(eList, indirect, gList) :

	# Collect unique gene names & make dict
	gList.sort()
	gDict = pp.createGeneMapping(gList)

	# Collect unique term names & make dict
#	eModList = eList[i for i in range(len(eList)) if elist[i,2] in indirect]
	eIdx = list()
	for i in range(len(eList)) :
		if eList[i,3] in indirect :
			eIdx.append(i)
	#end if
	eModList = eList[eIdx,:]
	tList = np.unique(eModList[:,0])
	tDict = pp.createGeneMapping(tList)

	# Allocate the gene-term matrix
	numG = len(gDict.keys())
	numT = len(tDict.keys())
	featMatrix = np.zeros( (numG, numT), dtype=featDType )

	# Track number of genes in each term
	tCount = np.zeros(numT, dtype=np.int32)

	# Place the weights into the matrix
	for edge in eModList :
		r = gDict[edge[1]]
		c = tDict[edge[0]]
		featMatrix[r,c] = edge[2]

		tCount[c] += 1
	#end loop

	# Save the feature matrix to a file
#	pp.setParamSaveTextCopy(True)
	pp.saveMatrixNumpy(featMatrix, 'featTerm_Orig', eDir, False)

	# Save the term names & counts to a file
	with open(eDir + 'featTerm_Names.txt', 'w') as fout :
		fout.write('{}\t{}'.format( tList[0], tCount[0] ))
		for i in range(1, len(tList)) :
			fout.write('\n{}\t{}'.format( tList[i], tCount[i]))
	#end with

	# Draw the term count distribution
	plt.hist(tCount, log=True, bins=200)
	plt.title('Distribution of Term Sizes')
	plt.xlabel('size of Term (# of Genes)')
	plt.ylabel('number of Terms')
	plt.savefig(eDir + 'featTerm_Dist.png')
#	plt.show()
	plt.close()

#end def ############### 




####### ####### ####### ####### 
# STANDALONE PROCEDURE

tstart = time.time()
print("\nCreating neighborhood features for {}".format(eName))

eDir = concatenatePaths(ePath, eName)

# GIVEN VARIABLES, if called from within preProcessing
geneHuman, keepGenes, loseGenes, keepEdges, indirEdges, thresh = pp.readKeepFile(eDir + 'keep.txt')
del geneHuman, loseGenes, keepEdges, thresh
edgeArray, nodeDict = pp.readEdgeFile(eDir + 'network.txt')
nodeDict, geneList = pp.createNodeLists(edgeArray, keepGenes)
del nodeDict

pp.setParamSaveTextCopy(saveAsText)

# Main Function call
print("Calling function createFeatureTerms() ...")
createFeatureTerms(edgeArray, indirEdges, geneList)


print("\nDone.\n")