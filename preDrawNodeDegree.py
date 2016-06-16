# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# Pre-Processing of the network
#	From the node-bin stats, draw degree distributions
#
# This assumes the pre-processing has already been run.
#	This script reads in node-degree.txt, and draws the
#	distribution of the node degree values for each of
#	edge types.
# ---------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import time



####### ####### ####### ####### 
# PARAMETERS

showFigures = True

# The network to use and directory path
#eName = 'all_v3beta_g2e9t0'
#ePath = '../Dropbox/mp/networks/'
eName = 'fakeNtwk00_g2e3t10'
ePath = 'networks/'

# # The number of bins to draw
numBinsMin = 4
numBinsMax = 50
#http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.hist

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()


# 1) Read in the file
if not eName.endswith('/') :
	eName = eName + '/'
if not ePath.endswith('/') :
	ePath = ePath + '/'
fName = ePath + eName + 'node-degree.txt'
print("Reading {} ...".format(fName))

#typeNames = list()
#nodeDegrees = np.zeros( ())

# get size of matrix
with open(fName, 'r') as fin :
	rows = -1
	for line in fin :
		rows += 1
	#end loop
	cols = len( line.split('\t') ) - 1
#end with

# read file ...
typeNames = list()
#geneNames = list()
nodeDegrees = np.zeros( (rows,cols) )
with open(fName, 'r') as fin :

	# Get the list of edge types
	header = fin.readline()
	header = header.rstrip()
	hv = header.split('\t')
	typeNames = hv[1:len(hv)]

	# Get the node-degree values
	r = -1
	for line in fin :
		r += 1
		line = line.rstrip()
		lv = line.split('\t')
		c = -1
		for value in lv[1:len(lv)] :
			c += 1
			nodeDegrees[r,c] = int(value)
#end with

#print(typeNames)
#print(nodeDegrees)



# 2) Determine the number of bins to use
if nodeDegrees.shape[0] > 10000 :
	useBins = numBinsMax
elif nodeDegrees.shape[0] > 50 :
	useBins = 10
else :
	useBins = numBinsMin
#end if



# 3) Draw the degree distribution images (and save)

for col in range(len(typeNames)) :
#for col in range(1) :

	thisName = typeNames[col]
	thisVals = nodeDegrees[:,col]

	plt.hist(thisVals, bins=useBins)
	plt.title('Node Degrees for type: {}'.format(thisName))
	plt.xlabel('Degree')
	plt.ylabel('Number of Nodes')

	plt.savefig(ePath + eName + 'nodeDeg-{}.png'.format(thisName))
	if showFigures :
		plt.show()
	plt.close()

#end loop



print("\nDone.\n")