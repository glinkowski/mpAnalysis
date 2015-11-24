# ----------------------------------------------------
# EXPERIMENT 16-2
# Goal: Find metapaths between two sets of genes
#
# From a sample file, select N source genes & N target
#   genes. Set a maxPathLength. For each source gene,
#   find all paths of maxPathLength or less. Take note
#   of what edge types are traversed, and how many
#   of that MetaPath exist.
# ----------------------------------------------------

import gfuncs as gf
import numpy as np
import random
from datetime import datetime



######## ######## ######## ########
# PARAMETERS

maxPathLength = 3
numSource = 1
numTarget = 1

outpath = "metapaths/"
outfile = "mp03-2.txt"

edgepath = "edgefiles/"
#edgefile = "h_gth_t0.edge.txt"
#edgefile = "fake.edge.txt"
edgefile = "hsa_dghmw.edge.txt"
delim = "\t"

samppath = "samples/"
sampname = "ACEVEDO_FGFR1_TARGETS_IN_PROSTATE_CANCER_MODEL"

######## ######## ######## ########



######## ######## ######## ########
# FUNCTION: DFS search for metapaths
# Returns: two dictionaries
#   mpathDict key = metapath (GO_term-prot_homol-GO_term)
#           value = number of occurrences of that path
#   lengthDict key = length of path (ie: 5 [edges])
#           value = number of paths of that length
######## ######## ######## ########
# FUNCTION: recursive helper for metapaths search
def findMetaPathsHelper(source, target, maxLength, network, mpathDict, lengthDict, edgeStack, pathStack) :

#    findMetaPathsHelper(source, target, maxLength, network, mpathDict, lengthDict, nodeQueue, currIdx)

#TODO: build this function!!

#    print edgeStack, pathStack
    # Pop (edge, node) from main edge stack
    mainPop = edgeStack.pop()

    # Add (edge, node) to current path stack
#    pathStack.append(mainPop)
    pathStack.append(mainPop[0])

    # Ending Conditions ########
    #   1) reached target node
    if (mainPop[1] == target) :
        # Record the path info
#        print "  reached target"
        # Update the metapaths dict
        mpkey = str()
        for i in pathStack :
#            print i
#            mpkey = mpkey + network[i[0],3] + "-"
            mpkey = mpkey + network[i,3] + "-"
        # remove last dash
        mpkey = mpkey.rstrip("-")
        # add to mpath dict
        if (mpkey in mpathDict.keys()) :
            mpathDict[mpkey] += 1
        else :
            mpathDict[mpkey] = 1
        #end loop

        # Update the length dict
        lengthDict[len(pathStack)] += 1

        # Remove newest edge from current path and return
        pathStack.pop()
##        return mpathDict, lengthDict
#        print "returning {}".format(edgeStack)
        return

    # 2) at maximum path length
    elif ( len(pathStack) >= maxLength) :
#        print "  reached max length"
        # Remove newest edge from current path and return
        pathStack.pop()
#        print "returning {}".format(edgeStack)
        return

    # 3) don't traverse same edge multiple times
#    elif (mainPop[0] in pathStack) :
    elif ( pathStack.count(mainPop[0]) > 1 ) :
        pathStack.pop()
        return
    #end if ########


    # Find all edges adjacent to current node
    cNode = mainPop[1]
    idxList = list()
    idxList.extend( np.where(network[:,0] == cNode)[0] )
    idxList.extend( np.where(network[:,1] == cNode)[0] )
    idxList = np.unique(idxList)
    idxList.sort()


    # Add new (edge, node) to the main edge stack
    numNew = 0
    for n in idxList :
        # skip the current edge
        if (n == mainPop[0]) :
            continue
        #end if
        node0 = network[n,0]
        node1 = network[n,1]
        numNew += 1
        if (node0 == mainPop[1]) :
            edgeStack.append( (n,node1) )
        else :
            edgeStack.append( (n,node0) )
    #end loop

    # Call self once for each added node
    for x in range(numNew) :
#        print "sending edges {}".format(edgeStack)
        findMetaPathsHelper(source, target, maxLength, network, mpathDict, lengthDict, edgeStack, pathStack)
    #end loop

    # Remove current (edge, node) and return
    pathStack.pop()
#    edgeStack.pop()
#    print "returning {}".format(edgeStack)
    return

#end def ######## ######## ########
######## ######## ######## ########
# FUNCTION: DFS search for metapaths
# Returns: two dictionaries
#   mpathDict key = metapath (GO_term-prot_homol-GO_term)
#           value = number of occurrences of that path
#   lengthDict key = length of path (ie: 5 [edges])
#           value = number of paths of that length
def findMetaPaths(source, target, maxLength, network) :

    # Initialize dicts to store path info
    mpathDict = dict()
    lengthDict = dict()
    for i in range(0,maxLength+1) :
        lengthDict[i] = 0
    #end loop

#    nodeStack = list()
#    nodeStack.append(source)
##TODO: add source to nodeQueue

    # main stack will contain (edge, node)
    edgeStack = list()
    # path stack will contain info on current path
    pathStack = list()


    # Find all edges adjacent to current node
#    cNode = mainPop[1]
    idxList = list()
    idxList.extend( np.where(network[:,0] == source)[0] )
    idxList.extend( np.where(network[:,1] == source)[0] )
    idxList = np.unique(idxList)
    idxList.sort()


    # Add each (edge, node) to the main edge stack
    for n in idxList :
        node0 = network[n,0]
        node1 = network[n,1]
        if (node0 == source) :
            edgeStack.append( (n,node1) )
        else :
            edgeStack.append( (n,node0) )
    #end loop

    # Call self till stack empty
    while edgeStack :
        findMetaPathsHelper(source, target, maxLength, network, mpathDict, lengthDict, edgeStack, pathStack)
    # end loop


    return mpathDict, lengthDict
#end def ######## ######## ########



######## ######## ######## ########
# BEGIN MAIN FUNCTION

print ""

# Get the Edges & Genes from the network file
print "Reading in the graph ..."
Edges, Nodes = gf.readGraphData(edgepath + edgefile, delim)

#TODO: write new gene function that will return all species
GeneSet = gf.readOnlyGenes(edgepath + edgefile, delim)

#print GeneSet

print "Extracting source & target genes from {}".format(sampname)
# Get source & target genes from the cancer sample
fnUp = sampname + "_UP.txt"
fnDn = sampname + "_DN.txt"
sampList = list()
fin = open(samppath + fnUp, "rb")
for line in fin :
    sampList.append(line.rstrip('\n'))
fin.close()
fin = open(samppath + fnDn, "rb")
for line in fin :
    sampList.append(line.rstrip('\n'))
fin.close()
sampList.sort()

#print "ready"

# Select the random N genes from cancer set
sourceList = list()
targetList = list()
for n in range(0,numSource) :
    sourceList.append(random.choice(sampList))
for n in range(0, numTarget) :
    targetList.append(random.choice(sampList))
#end loop
#while (len(sourceList) < numSource) :
#    thisOne = random.choice(sampList)
##    print thisOne, len(sourceList)
#    if (thisOne in GeneSet) :
#        print thisOne, numSource, len(sourceList)
#        sourceList.append(thisOne)
#while (len(targetList) < numTarget) :
#    thisOne = random.choice(sampList)
#    if (thisOne in GeneSet) :
#        targetList.append(thisOne)
##end loop

#print sourceList, targetList

# Identify genes in sample that aren't in network
missing = set(sourceList).difference(GeneSet)
missing = missing.union( set(targetList).difference(GeneSet) )
#print missing



print "Saving metapaths to file: {}".format(outfile)


timestart = datetime.now()

fout = open(outpath + outfile, "wb")
fout.write("{}\n".format(sampname))
fout.write("{}\n".format(edgefile))
fout.write("maxPath: {}\n".format(maxPathLength))
fout.write("sources: {}\n{}\n".format(numSource, sourceList))
fout.write("targets: {}\n{}\n".format(numTarget, targetList))
fout.write("\n")
for s in sourceList :
    for t in targetList :

        fout.write("from {} to {}\n".format(s,t))
        fout.close()

        time1 = datetime.now()
        mpathDict, lengthDict = findMetaPaths(s, t, maxPathLength, Edges)
        time2 = datetime.now()

#        print mpathDict

        fout = open(outpath + outfile, "ab")
        keyList = mpathDict.keys()
        keyList.sort()
        if not keyList :
            fout.write(" zero paths \n")
        for key in keyList :
            fout.write(" {:>3}\t{}\n".format(mpathDict[key], key))
        #end loop
        fout.write("\n")

        print "{} to {} has {} paths; time = {}".format( s,t,sum(mpathDict.values()),(time2 - time1) )
    #end loop
#end loop

timestop = datetime.now()
fout.write("\n{}".format(timestop - timestart))
fout.close()



print "\nDone.\n"