# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# Find all the metapaths in the network
#
# This is a basic demo. Load a sample file and rank all the
#   paths found between the genes in that sample. Do this
#   by comparing against random samples of the same size to
#   see which paths are more pronounced in the sample.
# Outline:
#   
# ---------------------------------------------------------

import mpFindFuncs as ff
import time



####### ####### ####### ####### 
# PARAMETERS

# The network to use and directory path
ename = 'fakeNtwk00_g2e3t10'
epath = 'networks/'

# The sample to test and path
sname = 'Fake00_sample01'
spath = 'samplesFake/'

# Where to store the output
oname = 'find01-' + ename + "-" + sname
opath = 'outputFake/'

# How many random samples to examine
numRand = 13

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()


# open the sample file
# get the list of paths available
# for each path, store:
#   this sample count (return left-out genes)
#   rand: mean, std dev, z-score
# output to file


# 1) Read in the sample file

print "Finding metapaths in sample: {}".format(sname)
sampGenes = ff.readSampleFiles(spath + sname, True, True)
print sampGenes


# 2) Identify the paths available

print "Checking what paths are available ..."
mpDict = ff.readKeyFile(epath, ename)
#print mpDict
#
mpList = mpDict.keys()
mpList.sort()
#print mpList




# FOR REFERENCE: copied from mpFindPaths00
## path to the metapath matrices
#mpPath = epath + ename + "_MetaPaths/"
## file containing the genes
#gfile = "genes.txt"
## Create mapping from gene name to row/col index
#geneDict = ff.readGeneFile(mpPath + gfile)
##print geneDict
## output file
#ofile = opath + oname + ".txt"
