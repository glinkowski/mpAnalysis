


import preProcFuncs as pp


####### ####### ####### ####### 
# PARAMETERS

ename = 'fakeNtwk00'
epath = '../networks/'

kfile = ename + '.keep.txt'


#mfile = ename + '.metapaths.txt'
#gfile = ename + '.genes.txt'
#rfile = ename + '.meta.gene_pairs.txt'
#pfile = ename + '.meta.path_names.txt'
delim = '\t'

####### ####### ####### ####### 


# read in the keep file
#   if not there, ask to have it created

print "Reading in the network:", ename
print "    reading the keep file", kfile
keepGenes, keepEdges, indirEdges = pp.readKeepFile(
    epath+kfile)

#print keepGenes
#print keepEdges
#print indirEdges




# read in the network to a matrix

# read in the corrections file
#   if not there, skip corrections
#   also, alert user

# ?
# change indirect edges to direct
# normalize weights

# create matrices and network
