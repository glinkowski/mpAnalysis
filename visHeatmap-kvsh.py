# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# View heatmap of connections between Known and Hidden genes
#
# Create a visualization -- heatmap -- of connections
#	between genes used in a prediction. Separate into three
#	groups: the Known genes, the Hidden/Concealed genes, and
#	a random selection of genes outside the sample.
# ---------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import visLibrary as vl



####### ####### ####### ####### 
# PARAMETERS

# Paths to network files & sample prediction files
nPath = '../Dropbox/mp/networks/'
nFolder = 'toy2_p3gz'
pPath = '../Dropbox/mp/output/'
pFolder = 'pred01-CAMPS_CO-002'

# Number of random genes to select
numRand = 200


#TODO: Can get the nFolder from scores.txt
if nFolder.endswith('/') :
	nDir = nPath + nFolder
else :
	nDir = nPath + nFolder + '/'
if pFolder.endswith('/') :
	pDir = pPath + pFolder
else :
	pDir = pPath + pFolder + '/'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION
print ""


# Get list of known genes
kGenes = vl.readFileColumnAsString(pDir+'known.txt', 0, 0)
#print len(kGenes)#, kGenes

# Get list of known genes
hGenes = vl.readFileColumnAsString(pDir+'concealed.txt', 0, 0)
#print len(hGenes)#, hGenes

# Get selection of random genes from network
allGenes = vl.readFileColumnAsString(nDir+'genes.txt', 0, 0)
#print len(allGenes)
sGenes = set(kGenes).union( set(hGenes) )
#print len(sGenes)
rGenes = vl.randSelectWithExclude(allGenes, sGenes, numRand)
#print len(rGenes)






# Numpy structured array to put genes in order:
# gene name, order (1, 2, 3 == rand, hidd, known), avg path count
