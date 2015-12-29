


import preProcFuncs as pp


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
print keepGenes
print loseGenes
print keepEdges
print indirEdges


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
	keepEdges)
print edgeArray


# Question: which should I keep: allGenes,
#   keepGenes, or loseGenes? I only need two.
allGenes = keepGenes + loseGenes
print allGenes
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


# save edge list, node dict, genes?
outname = pp.createModEdgeFileName(ename, keepEdges,
	keepGenes, thresh)
pp.writeModEdgeFilePlus(epath, outname, nodeDict,
	geneList, edgeArray)





# create the primary matrices
#   change indirect edges to direct
matrixList, matrixNames = pp.createMatrixList(edgeArray,
	keepEdges, indirEdges, geneList, outname)


# skip? - save updated network & node dict
#   no current use for it
#   BUT ... might be useful for other things
#       like DFS verification

# ? save genes to file
# ?

# create matrices and network
# save path types to file



print "\nDone.\n"