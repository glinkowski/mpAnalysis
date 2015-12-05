# ----------------------------------------------------
# author: Greg Linkowski
# project: metapaths for KnowEnG
#
# DFS search for metapaths between nodes
#
# This is to test the functions built from the original
#   mpDFSWalk script.
# From a sample file, select N source genes & N target
#   genes. Set a maxPathLength. For each source gene,
#   find all paths of maxPathLength or less. Take note
#   of what edge types are traversed, and how many
#   of that MetaPath exist.
# ----------------------------------------------------

import mpfuncs as mp
import numpy as np
import random




######## ######## ######## ########
# PARAMETERS

maxPathLength = 3
numSource = numTarget = 1

outpath = "../output/"
outfile = "mpDFS-test02.txt"

edgename = "toy_hsa_c"
edgefile = "../networks/" + edgename + ".edge.txt"
genefile = "../networks/" + edgename + ".genes.txt"
delim = "\t"

#samppath = "../samples/"
#sampname = "ACEVEDO_FGFR1_TARGETS_IN_PROSTATE_CANCER_MODEL"
sampname = "Random_from_all_genes_in_network"

######## ######## ######## ########




######## ######## ######## ########
# BEGIN MAIN FUNCTION

print ""

# Get the Edges & Nodes from the network file
print "Reading in the graph ..."
Edges, Nodes = mp.readEdgeFile(edgefile, delim)
del Nodes # not needed



# Get the Genes known to be in the network
geneList = list()
gf = open(genefile, "rb")
for line in gf:
    lv = line.split(delim)
    geneList.append(lv[0])
#end loop

# Select N random genes from gene set
sourceList = list()
for n in range(0,numSource) :
    sourceList.append(random.choice(geneList))
#end loop
targetList = sourceList
del geneList # no longer needed



# Run the DFS metapath search & save output
print "Saving metapaths to file: {}".format(outfile)

mp.outputMetaPathsDFS(outpath+outfile, maxPathLength, sampname, 
    edgename, Edges, sourceList, targetList)



print "\nDone.\n"