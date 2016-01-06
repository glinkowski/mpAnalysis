


import preProcFuncs as pp

import numpy as np


####### ####### ####### ####### 
# PARAMETERS

ename = 'fakeNtwk00'
epath = '../networks/'

kfile = ename + '.keep.txt'
efile = ename + '.edge.txt'
cfile = ename + '.correct.txt'


#mfile = ename + '.metapaths.txt'
#gfile = ename + '.genes.txt'
#rfile = ename + '.meta.gene_pairs.txt'
#pfile = ename + '.meta.path_names.txt'
#delim = '\t'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION


####
# read in the keep file
#   if not there, ask to have it created
print "\nReading in the network:", ename
print "    reading the keep file", kfile
keepGenes, loseGenes, keepEdges, indirEdges, thresh = pp.readKeepFile(epath+kfile)
#print keepGenes
#print loseGenes
#print keepEdges
#print indirEdges


####
# read in the network to a matrix
print "    reading in the edge file", efile
edgeArray, nodeDict = pp.readEdgeFile(epath+efile)
#print edgeArray


####
# read in the corrections file
#   if not there, skip corrections
#   also, alert user
print "Applying corrections to spelling ..."
#edgeArray = pp.applyCorrections(edgeArray, epath+cfile)
pp.applyCorrections(edgeArray, epath+cfile)
#print edgeArray


# normalize weights
print "Normalizing weights by edge type ..."
pp.applyNormalization(edgeArray, 0)
#print edgeArray


# threshold according to edge weights
#thresh = 3
print "Thresholding weights at {}".format(thresh)
edgeArray = pp.applyThreshold(edgeArray, thresh)
#print edgeArray


# throw out specified genes, edges
edgeArray = pp.applyKeepLists(edgeArray, loseGenes,
	keepEdges, indirEdges)
#print edgeArray


# Question: which should I keep: allGenes,
#   keepGenes, or loseGenes? I only need two.
allGenes = keepGenes + loseGenes
#print allGenes
#allGenes = keepGenes
#for regex in loseGenes :
#    allGenes.append(regex)
##end loop
#print allGenes

# create an updated nodeDict from modified edge list
nodeDict, geneList = pp.createNodeLists(edgeArray, allGenes)
# create an updated keepGenes from modified edge list
#geneList = createGeneList(edgeArray, keepGenes)
#print nodeDict
#print geneList


#this = set()
#this.add('apple')
#this.add('orange')
#print this
#this.add('apple')
#print this
#that = set()
#that.add('banana')
#those = this.union(that)
#print those


# skip? - save updated network & node dict
#   no current use for it
#   BUT ... might be useful for other things
#       like DFS verification
# ? save genes to file
# ?

# save edge list, node dict, genes?
outname = pp.createModEdgeFileName(ename, keepEdges,
	keepGenes, thresh)
print "Saving modified network to {}.edge.txt".format(outname)
pp.writeModEdgeFilePlus(epath+"modified/", outname,
	nodeDict, geneList, edgeArray)


# create the primary matrices
#   change indirect edges to direct
print "Creating the primary gene-gene matrices ..."
matrixList, matrixNames = pp.createMatrixList(edgeArray,
	keepEdges, indirEdges, geneList, nodeDict)#, epath, outname)
print matrixNames
print matrixList[0]

# save the primary matrices
primpath = epath + outname + "_Primaries/"
pp.saveMatrixList(matrixList, matrixNames, geneList, primpath)


# create matrices (? and network)
# save path types to file
mpPath = epath + outname + "_MetaPaths/"
# Create the 1-step paths
#primNames = matrixNames #NOPE! This just makes a pointer
primNames = list()
mDict = dict()	# map indices to the paths
mNum = 0
for name in matrixNames :
	primNames.append(name)
	mDict[name] = mNum
	mNum += 1
#end loop
# Create the 2-step paths
for i in range(0, len(primNames)) :
	for j in range(i, len(primNames)) :
		newM = np.dot(matrixList[i], matrixList[j])
		name1 = primNames[i] + "-" + primNames[j]
		name2 = primNames[j] + "-" + primNames[i]
		if name1 == name2 :
			matrixList.append(newM)
			matrixNames.append(name1)
			mDict[name1] = mNum
		else :
			matrixList.append(newM)
			matrixNames.append(name1)
			mDict[name1] = mNum
			mDict[name2] = mNum
		#end if

#		newM = np.multiply(matrixList[i], matrixList[j])
#		name = primNames[i] + "-" + primNames[j]
#		matrixList.append(newM)
#		matrixNames.append(name)

		mNum += 1
	#end loop
#end loop

#pp.saveMatrixList(matrixList, matrixNames, geneList, mpPath)
pp.saveMatrixListPlus(matrixList, mDict, geneList, mpPath)

# Create the 3-step paths
# Create the 4-step paths

#pp.saveMatrixList(matrixList, matrixNames, geneList, mpPath)
#pp.saveMatrixListPlus(matrixList, mDict, geneList, mpPath)


print "\nDone.\n"