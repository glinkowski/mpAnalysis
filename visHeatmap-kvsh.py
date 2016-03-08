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
# NOTE: This only looks at the primary edge types, or
#	metapaths of length 1.
# ---------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import visLibrary as vl
import mpLibrary as mp
import preProcFuncs as pp



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
#eTypes = vl.readFileColumnAsString(nDir+'edges.txt', 0, 0)

# Get the primary path matrices & names
#mpDir = nDir[0:-1]+'-Primaries/'
pathDict = pp.readPrimaryMatrices(nPath, nFolder)
eTypes = pathDict.keys()
eTypes.sort()

# Define the matrix to hold the path counts (the heatmap)
hmap = np.zeros([len(gOrder), len(gOrder)])


# Identify the hottest genes within gOrder
indices = [geneDict[g] for g in gOrder['name']]
indices.sort()
#print indices


for et in eTypes :
#et = 'prot_homol'


	# reduce the path matrix to the desired columns
	#avgCounts = np.mean( pathDict[et][:,indices], axis=1 )
	sumCounts = np.sum( pathDict[et][:,indices], axis=1 )
	#print len(avgCounts), pathDict[et].shape[0]
	# for each gene in gOrder, get the avg path count w/in this set
	#gOrder['pcount'] = [avgCounts[geneDict[g]] for g in gOrder['name']]
	gOrder['pcount'] = [sumCounts[geneDict[g]] for g in gOrder['name']]

	gOrder.sort(order=['order', 'pcount'])
	#print gOrder

	# Fill the heatmap with values from path matrix
	for x in xrange(len(gOrder)) :
		for y in xrange(x,len(gOrder)) :

			hmap[x,y] = (pathDict[et][ geneDict[gOrder['name'][x]],
				geneDict[gOrder['name'][y]] ])

			if x != y :
				hmap[y,x] = hmap[x,y]
	#end loop
	hmax = np.amax(hmap)


#TODO: Should I apply a scaling to better visualize this? A log transform?
	# Apply log scaling
#	hmap = np.log(hmap + 1)


	# Plot the figure
	plt.pcolor(hmap, cmap=plt.cm.Reds)
	plt.suptitle('Heatmap of '+pFolder+' for edge '+et)
	plt.title('max count = {}'.format(int(hmax)))
	# omit whitespace & reverse the order of the axes
	plt.axis([0, len(gOrder), 0, len(gOrder)])
	plt.gca().invert_xaxis()
	plt.gca().invert_yaxis()
	# label the axis labels and ticks
#	plt.yticks(range(len(gOrder)), gOrder['name'], fontsize=5)
	plt.yticks([len(gOrder) - len(sGenes), len(gOrder) - len(kGenes),
		len(gOrder)], ['outside set', 'hidden', 'known'],
		ha='right', va='bottom', rotation='70')
	plt.ylabel('genes: ascending by avg connections within this set')
	plt.xticks([len(gOrder) - len(sGenes), len(gOrder) - len(kGenes),
		len(gOrder)], ['outside set', 'hidden', 'known'], ha='left', rotation='-20')
	plt.plot([0,len(gOrder)],[len(gOrder) - len(sGenes),len(gOrder) - len(sGenes)], 'lightgrey')
	plt.plot([0,len(gOrder)],[len(gOrder) - len(kGenes),len(gOrder) - len(kGenes)], 'lightgrey')
	plt.plot([len(gOrder) - len(sGenes),len(gOrder) - len(sGenes)],[0,len(gOrder)], 'lightgrey')
	plt.plot([len(gOrder) - len(kGenes),len(gOrder) - len(kGenes)],[0,len(gOrder)], 'lightgrey')
	plt.subplots_adjust(left=0.2)
	print "  saving Heatmap of {}".format(et)
	plt.savefig(pDir+'Heatmap_'+et+'.png')
#	plt.show()
	plt.close()

#end loop



# Get the top paths from this prediction
topPaths, topPathScores = vl.getTopRankedItems(pDir+'ranked_paths.txt', 10, 0)
#print topPaths, topPathScores

pathTuples = mp.readKeyFile(nPath, nFolder)

matrixSize = (pathDict[et]).shape[0]
#print matrixSize
for pt in topPaths :

	matrix = mp.getPathMatrix(pathTuples[pt], nPath, nFolder, matrixSize)

	# reduce the path matrix to the desired columns
	sumCounts = np.sum( matrix[:,indices], axis=1 )
	# for each gene in gOrder, get the avg path count w/in this set
	gOrder['pcount'] = [sumCounts[geneDict[g]] for g in gOrder['name']]

	gOrder.sort(order=['order', 'pcount'])

	# Fill the heatmap with values from path matrix
	for x in xrange(len(gOrder)) :
		for y in xrange(x,len(gOrder)) :

			hmap[x,y] = (matrix[ geneDict[gOrder['name'][x]],
				geneDict[gOrder['name'][y]] ])

			if x != y :
				hmap[y,x] = hmap[x,y]
	#end loop
	hmax = np.amax(hmap)


	# Plot the figure
	plt.pcolor(hmap, cmap=plt.cm.Reds)
	plt.suptitle('Heatmap of '+pFolder+' for '+pt)
	plt.title('max count = {}'.format(int(hmax)))
	# omit whitespace & reverse the order of the axes
	plt.axis([0, len(gOrder), 0, len(gOrder)])
	plt.gca().invert_xaxis()
	plt.gca().invert_yaxis()
	# label the axis labels and ticks
#	plt.yticks(range(len(gOrder)), gOrder['name'], fontsize=5)
	plt.yticks([len(gOrder) - len(sGenes), len(gOrder) - len(kGenes),
		len(gOrder)], ['outside set', 'hidden', 'known'],
		ha='right', va='bottom', rotation='70')
	plt.ylabel('genes: ascending by avg connections within this set')
	plt.xticks([len(gOrder) - len(sGenes), len(gOrder) - len(kGenes),
		len(gOrder)], ['outside set', 'hidden', 'known'], ha='left', rotation='-20')
	plt.plot([0,len(gOrder)],[len(gOrder) - len(sGenes),len(gOrder) - len(sGenes)], 'lightgrey')
	plt.plot([0,len(gOrder)],[len(gOrder) - len(kGenes),len(gOrder) - len(kGenes)], 'lightgrey')
	plt.plot([len(gOrder) - len(sGenes),len(gOrder) - len(sGenes)],[0,len(gOrder)], 'lightgrey')
	plt.plot([len(gOrder) - len(kGenes),len(gOrder) - len(kGenes)],[0,len(gOrder)], 'lightgrey')
	plt.subplots_adjust(left=0.2)
	print "  saving Heatmap of {}".format(pt)
	plt.savefig(pDir+'Heatmap_'+pt+'.png')
#	plt.show()
	plt.close()

#end loop



print "\nDone.\n"