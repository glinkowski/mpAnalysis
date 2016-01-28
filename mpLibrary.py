# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# functions used by the Metapath Analysis scripts
#
# ---------------------------------------------------------

import os.path
from os import listdir
import sys
import numpy as np
import random



######## ######## ######## ######## 
# PARAMETERS

# Seed the random functions, if desired
random.seed(42)

# File names to use when partitioning gene set:
fnKnownGenes = 'known.txt'
fnConcealedGenes = 'concealed.txt'
fnLeftOutGenes = 'ignored.txt'

## Data type used when loading edge file to memory:
#nodeDT = np.dtype('a30')

######## ######## ######## ######## 






######## ######## ######## ########
# Function: save a list to a text file
# Input:
#   path, str: path to save the file
#   name, str: name of file to save
#   theList, list of str - list of items to save
#       ASSUMPTION: list is already properly ordered
# Returns: nothing
# Creates: file containing ordered list of items
def saveListToText(path, name, theList) :

    # If folder doesn't exist, create it
    if not os.path.exists(path) :
        os.makedirs(path)
    #end if

    theFile = open(path + name, 'wb')
    firstline = True
    for item in theList :
        if firstline :
            firstline = False
        else :
            theFile.write("\n")
        #end if
        theFile.write("{}".format(item))
    #end if
    theFile.close()

    return
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Take the list of genes from the sample, return
#   a list of genes in the network, and left-out genes
# Input ----
#   path, str: path to the network files
#   name, str: name of the network to use
#   geneList, str list: genes in the sample
# Returns ----
#   inList & outList, str list: which genes from the
#       given list are found in network & which aren't
def checkGenesInNetwork(path, name, geneList) :

    fname = path + name + '_MetaPaths/genes.txt'

    # ERROR CHECK: verify file exists
    if not os.path.isfile(fname) :
        print ( "ERROR: Specified file doesn't exist:" +
            " {}".format(fname) )
        sys.exit()
    #end if

    # Read in the genes from the file
    geneSet = set()
    fg = open(fname, 'rb')
    for line in fg :
        line = line.rstrip()
        geneSet.add(line)
    #end loop
    fg.close()

    # The items to return
    inList = list()
    outList = list()

    # Sift through the sample
    for gene in geneList :
        if gene in geneSet :
            inList.append(gene)
        else :
            outList.append(gene)
        #end if
    #end loop

    inList.sort()
    outList.sort()
    return inList, outList
#end def ######## ######## ######## 



######## ######## ######## ########
# Function: Given a set of sample genes, 
# Pre-req: checkGenesInNetwork(path, name, geneList)
# Input ----
#   nPath, str: path to the network files
#   nName, str: name of the network to use
#   oPath, str: directory to create output
#   sample, str list: all genes in the sample
#   percent, int/float: percent of sample to conceal
# Returns ----
#   kGenes, str list: genes to be used as "known" test set
#   hGenes, str list: genes to be used as concealed set
# Creates ----
#   outPath/known.txt: list of known genes
#   outPath/concealed.txt: list of hidden genes
#   outPath/ignored.txt: list of genes not in network
def partitionSample(nPath, nName, oPath, sample, percent) :

    # Convert percent to fraction, if needed
    if (percent > 0) and (percent < 1) :
        # Assume percent is given as float (0, 1)
        pHide = percent
    elif (percent < 100) :
        # Assume percent is given as int [1, 100)
        pHide = percent / float(100)
    else :
        print ("ERROR: Invalid value given as the percent" +
            " of sample to hide.")
    #end if

    # Remove any genes not in the network
    inGenes, outGenes = checkGenesInNetwork(nPath, nName, sample)

    # Number of genes to keep vs hide
    numHide = int(len(inGenes) * pHide)
    numKeep = len(inGenes) - numHide

    # Separate the known and hidden genes
    random.shuffle(inGenes)
    kGenes = inGenes[0:numKeep]
    hGenes = inGenes[numKeep:len(inGenes)]

    # Save the gene lists to files
    kGenes.sort()
    saveListToText(oPath, fnKnownGenes, kGenes)
    hGenes.sort()
    saveListToText(oPath, fnConcealedGenes, hGenes)
    saveListToText(oPath, fnLeftOutGenes, outGenes)

    return kGenes, hGenes
#end def ######## ######## ########



