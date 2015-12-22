


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


####
# read in the keep file
#   if not there, ask to have it created

print "\nReading in the network:", ename
print "    reading the keep file", kfile
keepGenes, keepEdges, indirEdges = pp.readKeepFile(
    epath+kfile)

#print keepGenes
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
thresh = 3
print "Thresholding weights at {}".format(thresh)
edgeArray = pp.applyThreshold(edgeArray, thresh)
print edgeArray

# skip? - save updated network & node dict
#   no current use for it

# ? save genes to file
# ?
# change indirect edges to direct

# create matrices and network
# save path types to file



print "\nDone.\n"