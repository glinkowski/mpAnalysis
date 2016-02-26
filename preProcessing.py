# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Pre-Processing of the network
#
# To speed the calculation of metapaths, this script
#	applies some pre-processing steps to the network.
# The end result is a set of 2-D matrices, each indicating
#	the number of paths of a certain type from one gene to
#	another, where each file corresponds to a specific
#	metapath.
# Outline:
#	1) Load the given network and fix known typos, extract
#		desired nodes & edges, and threshold edge weights
#	2) Save the modified network (w/ node-index dict)
#	3) Find the primary (single-step) path matrices
#	4) Use those to calculate the rest of the metapaths
# ---------------------------------------------------------

import preProcFuncs as pp
import time



####### ####### ####### ####### 
# PARAMETERS

# The network to use and directory path
useRealData = True

if not useRealData :
#	ename = 'fakeNtwk00_g2e3t10'
	ename = 'fakeNtwk01'
	epath = 'networks/'
else :
#	ename = all_v1
	ename = 'toy2_hsa'
	epath = '../Dropbox/mp/networks/'
#end if


# Maximum number of steps in the calculated metapaths
mpDepth = 4

kfile = ename + '.keep.txt'
efile = ename + '.edge.txt'
cfile = ename + '.correct.txt'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()



# 1a) Read in the network from file(s)

# Read in the keep file
#	if not there, ask to have it created
print "\nReading in the network:", ename
print "    reading the keep file", kfile
geneHuman, keepGenes, loseGenes, keepEdges, indirEdges, thresh = pp.readKeepFile(epath+kfile)
#print keepGenes
#print keepEdges
#print indirEdges

# Read in the network (edge file) to a matrix
print "    reading in the edge file", efile
edgeArray, nodeDict = pp.readEdgeFile(epath+efile)
print "    initial edgeArray size: {} bytes".format(edgeArray.nbytes)


# Read in the corrections file
#   if not there, skip corrections, alert user
print "Applying corrections to spelling ..."
pp.applyCorrections(edgeArray, epath+cfile)



# 1b) Apply normalization, thresholding, discard unwanteds

# Normalize weights
print "Normalizing weights by edge type ..."
pp.applyNormalization(edgeArray, 0)
# Apply threshold to normalized edges
print "Thresholding weights at {}".format(thresh)
edgeArray = pp.applyThreshold(edgeArray, thresh)

# Discard specified genes, edges
edgeArray = pp.applyKeepLists(edgeArray, loseGenes,
	keepEdges, indirEdges)
print "    final edgeArray size: {} bytes".format(edgeArray.nbytes)


print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 2) Save the modified network

# Create an updated nodeDict from modified edge list
nodeDict, geneList = pp.createNodeLists(edgeArray, keepGenes)

# Save the edge list, node dict, gene list
outname = pp.createModEdgeFileName(ename, keepEdges,
	keepGenes, thresh)
print "Saving modified network to {}.network.txt".format(outname)
pp.writeModEdgeFilePlus(epath, outname,
	nodeDict, geneList, edgeArray)

# Save the node-binning stats for the network
pp.saveSelectGeneDegrees(epath, outname, edgeArray, geneList, geneHuman)


print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 3) Find the primary path matrices

# Change indirect edges to direct
print "Creating the primary gene-gene matrices ..."
matrixList = pp.createMatrixListNoBinning(edgeArray,
	keepEdges, indirEdges, geneList, nodeDict)

# Save the primary matrices
primpath = epath + outname + "_Primaries/"
print "    ... saving to: {}".format(primpath)
pp.saveMatrixList(matrixList, geneList, primpath)

print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



# 4) Calculate the specified metapaths

print "Determining metapaths, up to length {}".format(mpDepth)
# folder/path to save output
mpPath = epath + outname + "_MetaPaths/"
print "    ... and saving to: {}".format(mpPath)
# Find paths up to length = mDepth
pp.createMetaPaths(matrixList, geneList, mpDepth, mpPath)

print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)



print "\nDone.\n"