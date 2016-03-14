# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# Draw area under the curve(s) for predictions
#
# Draw the area under the ROC and the PR curve for a given
#	prediction attempt.
# ---------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import visLibrary as vl
import mpLibrary as mp
import preProcFuncs as pp



####### ####### ####### ####### 
# PARAMETERS

# Paths to network files & sample prediction files
pPath = '../Dropbox/mp/output/'
pFolder = 'pred01-CAMPS_CO-006'


if pFolder.endswith('/') :
	pDir = pPath + pFolder
else :
	pDir = pPath + pFolder + '/'

# Data type used by preProcessing when reading in nodes
nodeDT = np.dtype('a30')


textDelim = '\t'
####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION
print ""


## Read in the hidden genes
#gHidden = vl.readFileColumnAsString(pDir+'concealed.txt', 1, 0)
##print gHidden

#
## Declare the confusion matrix
#rows, colMin, colMax = vl.countLinesInFile(pDir+'ranked_genes.txt')
##print rows, colMin, colMax
#posActual = len(gHidden)
#negActual = rows - posActual
#confusion = np.zeros([2,rows])	# TP, FP
#
## Read in the ranked genes
#gHidSet = set(gHidden)
#
#fin = open(pDir+'ranked_genes.txt')
#col = 0
#for line in fin :
#
#	line = line.rstrip()
#	lv = line.split(textDelim)
#
#	if lv[1] in gHidSet :
#		confusion[0,col] += 1
#	else :
#		confusion[1,col] += 1
#	#end if
#
#	col += 1
##end loop
#fin.close()

#
## Convert into Recall (FPR), TPR, & Precision
## TPR = TP / allP = confusion[0,i] / posActual
## FPR = FP / allN = confusion[1,i] / negActual
## recall = TPR
## precision = TP / TP + FP = confusion[0,i] / (confusion[0,i] + confusion[1,i])
#recall = confusion[0,i] / posActual		# aka TPR
#FPR = confusion[1,i] / negActual
#precision = confusion[0,i] / (confusion[0,i] + confusion[1,i])


# Get the curve statistics from the ranked_genes file
FPR, recall, precision = vl.getAUCstats(pDir)
print recall
print FPR
print precision

# plot the results
