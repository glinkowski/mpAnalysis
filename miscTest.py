# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# This is a test script to call and verify different functions
#
# ---------------------------------------------------------

import mpFindFuncs as ff
#import preProcFuncs as pp
import numpy as np
import re
import os


humanRegex = ['ENSG', 'LGI_']
keep = list()
for gene in ['ENSG09483', 'LGI_43287', 'tomato'] :
	# list of matches to be found
	ma = list()
	for exp in humanRegex :
		ma.append( re.match(exp, gene) )
	if any(match != None for match in ma) :
		keep.append(gene)
#end loop
print keep

#keep = list()
#for gene in ['ENSG09483', 'LGI_43287', 'tomato'] :
#	if any(match != None for match in re.match(humanRegex, gene)) :
#		keep.append(gene)
##end loop
#print keep


nPath = 'networks/'
nName = 'fakeNtwk01'
fname = nPath+nName+'/node-degree.txt'
# ERROR CHECK: rename if file exists
if os.path.isfile(fname) :
	print ( "WARNING: Specified file already exists:" +
		" {}".format(fname) )
	i = 0
	while os.path.isfile(fname) :
		fname = nPath+nName+"/node-degree{:03d}.txt".format(i)
		i += 1
	print "    using new file name: {}".format(fname)
#end if