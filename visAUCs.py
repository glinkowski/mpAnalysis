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
pFolder = 'pred01-CAMPS_CO-002'


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

