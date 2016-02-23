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
import mpLibrary as mp



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

# Data type used by preProcessing when reading in nodes
nodeDT = np.dtype('a30')

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


# Define the array used to sort the genes
# Numpy structured array to put genes in order:
# gene name, order (1, 2, 3 == rand, hidd, known), avg path count

gOrder = np.recarray( (len(rGenes) + len(sGenes)),
	dtype=[('name', nodeDT), ('order', 'i4'), ('pcount', 'f4')] )

row = 0
for g in rGenes :
	gOrder[row] = (g, 1, 0)
	row += 1
for g in hGenes :
	gOrder[row] = (g, 2, 0)
	row += 1
for g in kGenes :
	gOrder[row] = (g, 3, 0)
	row += 1
#end if

geneDict = mp.readFileAsIndexDict(nDir+'genes.txt')
eTypes = vl.readFileColumnAsString(nDir+'edges.txt', 0, 0)

mpDir = nDir[0:-1]+'-Primaries/'
#TODO: create mapping of etypes to matrix files


#TODO: define path count array

#for et in eTypes :

#TODO: create array of path counts, sort gOrder by order+pcount
#TODO: output a heatmap image to pDir





