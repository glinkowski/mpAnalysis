# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# Find all the metapaths in the network
#
# This is just a preliminary verification step. Within a
#   given network, find all the paths of each type that
#   occur between all the genes. Output this data to file.
# Outline:
#   
# ---------------------------------------------------------

import mpFindFuncs as ff
import time



####### ####### ####### ####### 
# PARAMETERS

# The network to use and directory path
ename = 'fakeNtwk00_g2t3h10'
epath = 'networks/'

# Where to store the output
oname = "test00_" + ename
opath = 'output/'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()


# get the list of genes in the network
# find & output all the metapaths among those genes

# path to the metapath matrices
mpPath = epath + ename + "_MetaPaths/"

# file containing the genes
gfile = "genes.txt"

# Create mapping from gene name to row/col index
geneDict = ff.readFileAsIndex(mpPath + gfile)
#print geneDict

# output file
ofile = opath + oname + ".txt"
