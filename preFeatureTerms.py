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
#import matplotlib.pyplot as plt
import numpy as np
import gzip
import time
import mpLibrary as mp




####### ####### ####### ####### 
# PARAMETERS

# The network to use and directory path
useFakeData = True
if useFakeData :
	eName = 'fakeNtwk00_g2e3t10'
	ePath = 'networks/'
else :
	eName = 'all_v3beta_g2e9t0'
	ePath = '../Dropbox/mp/networks/'
#end if

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

	#1 From the keep file, get the gene regex
	#	alt: get the indirect types
	#2 read in the edge file, skip direct edges
	#  read in the gene index dict

	#3 collect unique term names
	#  create gene-term matrix
	#4 For each line in the edge list
	#	track total genes in term
	#	place weights into matrix
	#5 output term list w/ total counts
	#  output matrix to file


	# Collect unique gene names & make dict
	gList.sort()
	gDict = pp.createGeneMapping(gList)

	# Collect unique term names & make dict
	eModList = eList[i for i in eList if eList[i,2] in indirect]
	tList = np.unique(eModList[:,2])
	tDict = pp.createGeneMapping(tList)

	# Allocate the gene-term matrix
	numG = len(gDict.keys())
	numT = len(tDict.keys())
	featMatrix = np.zeros( (numG, numT), dtype=featDType )

	# Track number of genes in each term
	tCount = np.zeros(numT, dtype=np.int32)

	# Place the weights into the matrix
	for edge in eList :
		# Skip if edge type is direct gene-gene
		if edge[3] not in indirect :
			continue

		r = gDict[edge[1]]
		c = tDict[edge[0]]
		featMatrix[r,c] = edge[2]

		tCount[c] += 1
	#end loop

	# Save the feature matrix to a file
	saveMatrixNumpy(featMatrix, 'featTerm_Orig', eDir, False)

	# Save the term names & counts to a file
	with open(eDir + 'featNeighbor_Names.txt', 'w') as fout :
		fout.write('{}\t{}'.format( tList[0], tCount[0] ))
		for i in range(1, len(tList)) :
			fout.write('\n{}\t{}'.format( tList[i], tCount[i]))
	#end with

#end def ############### 




####### ####### ####### ####### 
# STANDALONE PROCEDURE

tstart = time.time()
print("\nCreating neighborhood features for {}".format(eName))

eDir = concatenatePaths(ePath, eName)

# GIVEN VARIABLES, if called from within preProcessing
geneHuman, keepGenes, loseGenes, keepEdges, indirEdges, thresh = pp.readKeepFile(eDir + 'keep.txt')
edgeArray, nodeDict = pp.readEdgeFile(eDir + 'network.txt')
nodeDict, geneList = pp.createNodeLists(edgeArray, keepGenes)
#TODO: delete what I don't need

# Main Function call
createFeatureTerms(edgeArray, indirEdges, geneList)


print("\nDone.\n")