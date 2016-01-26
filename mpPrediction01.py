# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Predict membership of gene sets using metapaths
#
# This code is to test the idea of using metapath analysis
#	as a one-class classifier. In this script, a portion
#	of a known set will be concealed, and the script will
#	attempt to predict the concealed members.
# Outline:
#	1) 
# ---------------------------------------------------------

import mpFindFuncs as ff
#import preProcFuncs as pp
import time



####### ####### ####### ####### 
# PARAMETERS

# The percent of genes to conceal
percHide = .25
# The number of predicted genes to return
numTopK = 100

# The network & sample sets to use (0 means fake)
useNtwk = 0

if useNtwk == 0 :
	ename = 'fakeNtwk00'
	epath = 'networks/'
	sname = 'sample01'
	spath = 'samplesFake/'
else :
	ename = 'toy2_hsa'
	epath = '../networks/'
	sname = 'XXX'
	spath = '../samples/'
#end if

# The path & file name to save output (& temp)
#oname = ''
opath = 'pred01-' + sname[0:7]

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()


# 1) Load a sample and hide some genes

# Read in the full sample
print "Analyzing metapaths in sample: {}".format(sname)
fullSample = ff.readSampleFiles(spath + sname, True, True)
# Remove any genes not in the network
inGenes, outGenes = ff.checkGenesInNetwork(epath,
	ename, fullSample)
print ("Of the {} sample genes,".format(len(sampGenes)) +
	" {} are in the network.".format(len(inGenes)) )

# select genes to test (function: opath, inGenes, outGenes, percHide)
# save to files: known, concealed, ignored


