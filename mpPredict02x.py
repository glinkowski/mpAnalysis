# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Approach 2: learn top paths from PathSim sum, find genes
#	Version 2, Step 3
#	Compile the Data
#
# The first approach tries to first find the most
#	important metapaths, then rank genes by similarity
#	along those. This approach instead finds the set
#	similarity along each metapath, then tries to learn
#	the most important paths.
#
# Step 3: After creating results from all of the attempts
#	to rank genes (regressions, etc), collect the results.
#	Create matrices to pull into Excel:
#	- sample-by-method, AUC values (ROC & Prec-Recall)
#	- sample-by-metapath, avg path weight after normalizing
# ---------------------------------------------------------

import mpLibrary as mp
#import time
import numpy as np
#import gzip
import visLibrary as vl
import matplotlib.pyplot as plt
import os
import sys



####### ####### ####### ####### 
# PARAMETERS

# folder containing the samples & results
#dDir = 'pred04-set02'
dDir = 'pred02a-batch-002'
dRoot = '../Dropbox/mp/output/'

# Network name & location
eName = 'all_v3beta_g2e9t0'
ePath = '../Dropbox/mp/networks/'


# verbose feedback ?
newVerbose = True

textDelim = '\t'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

#tstart = time.time()
print("")

mp.setParamVerbose(newVerbose)



# 1) Load the path names

# Read the network location from parameters file
if dDir.endswith('/') :
	dPath = dRoot + dDir
else :
	dPath = dRoot + dDir + '/'
#end if

print("Creating the gene-index dictionary.")
geneDict = mp.readGenesFile(ePath, eName)
geneList = list(geneDict.keys())
geneList.sort()

print("Reading in the path names.")
pathDict = mp.readKeyFile(ePath, eName)
pathList = mp.removeInvertedPaths(pathDict)
del pathDict

# Make new pathDict to give index of pathList
pathDict = dict()
idx = -1
for item in pathList :
	idx += 1
	pathDict[item] = idx
#end loop



# 2) Get the list of samples
dSubDirs = mp.getSubDirectoryList(dRoot+dDir)
sampList = list()
for sd in dSubDirs :
	sv = sd.split('/')
	sampList.append(sv[-2])
#end loop



# 3) Create the matrices to hold the results

# Get list of files in the subdirectory
fileList = os.listdir(dSubDirs[0])
fileList.sort()

# ERROR CHECK: ensure ranked_genes & ranked_paths files exist
rg = 0
rp = 0
for item in fileList :
	iv = item.split('/')
	if iv[-1].startswith('ranked_genes') :
		rg += 1
	elif iv[-1].startswith('ranked_paths') :
		rp += 1
#end loop
if rg < 1 :
	print("ERROR: No ranked_genes... files exist in {}".format(dSubDirs[0]))
	sys.exit
#end if

# Get list of the methods used from the file names
methodList = list()
for item in fileList :
	iv = item.split('/')
	fn = iv[-1]
	if fn.startswith('ranked_genes') :
		fn = fn[0:-4]
		fv = fn.split('-')
		if (fv[1] != 'votingAll') :
			methodList.append(fv[1] + '-' + fv[2])
#end loop
methodList.sort()


resultsROC = np.zeros( (len(methodList), len(dSubDirs)), dtype=np.float64)
resultsPR = np.zeros( (len(methodList), len(dSubDirs)), dtype=np.float64)
#resultsPaths = np.zeros( (len(pathList), len(dSubDirs)), dtype=np.float64)
#resultsGenes = np.zeros( (len(geneList), len(dSubDirs)), dtype=np.float64)



# 4) Create the Area Under the Curve tables
print("Finding AUCs for each sample ...")

col = -1
#for sd in dSubDirs[0:1] :
for sd in dSubDirs :
	col += 1

	# Get data relating to each method
	row = -1
	for m in methodList :
		row += 1

		fn = 'ranked_genes-' + m + '.txt'
		FPR, recall, precision, numHid = vl.getAUCstats(sd, fn)

		# Calculate (approximate) are under the ROC curve
		areaROC = 0
		for r in recall :
			areaROC += (r / len(recall))
		#end loop

		# Calculate (approximate) are under the PR curve
		areaPR = 0
		for p in precision :
			areaPR += (p / len(precision))
		#end loop

		# save data into the matrix
		resultsROC[row,col] = areaROC
		resultsPR[row,col] = areaPR


#TODO: this should be a func
		# Save the AUC figure(s)
		outName = sd + 'AUC-' + m + '.png'
		fig = plt.figure()

		# Plot the ROC curve
		plt.subplot(1, 2, 1)
		plt.plot(FPR, recall)
		plt.plot([0,1], [0,1], 'lightgrey')
		plt.xlabel('False Positive Rate')
		plt.ylabel('True Positive Rate')
		plt.axis([0, 1, 0, 1])

		# Plot the Precision-Recall curve
		plt.subplot(1, 2, 2)
		plt.plot(recall, precision)
		plt.xlabel('Recall')
		plt.ylabel('Precision')
		plt.axis([0, 1, 0, 1])

		# Final touches
		sv = sd.split('/')
		sn = sv[-2]
		plt.suptitle(sn + '\n{}, concealed = {}'.format(m, numHid) +
			', ROC area = {:.3}'.format(areaROC))
		plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.4, hspace=None)

		# Save the figure
		plt.savefig(outName)
		plt.close()
#end loop



# 5) Output to file(s)
print("Writing the AUC tables (ROC & PR) ...")

#TODO: write this into a func

# Write the AUC tables to file(s)
with open(dPath + 'results-AUC_ROC.txt', 'w') as fout :
	fout.write('Area Under ROC Curve, per-sample')
	fout.write('\nnetwork:{}{}'.format(textDelim, eName))
#	fout.write('\nfolds:{}{}'.format(textDelim, numFolds))
	fout.write('\n\n')

	for j in range(len(sampList)) :
		fout.write('{}{}'.format( sampList[j], textDelim ))
	#fout.write('\n')

	for i in range(len(methodList)) :
		fout.write('\n')
		for j in range(len(sampList)) :
			fout.write('{}{}'.format( resultsROC[i,j], textDelim ))
		fout.write('{}'.format(methodList[i]))
	#end loop
#end with

with open(dPath + 'results-AUC_PR.txt', 'w') as fout :
	fout.write('Area Under PR Curve, per-sample')
	fout.write('\nnetwork:{}{}'.format(textDelim, eName))
#	fout.write('\nfolds:{}{}'.format(textDelim, numFolds))
	fout.write('\n\n')

	for j in range(len(sampList)) :
		fout.write('{}{}'.format( sampList[j], textDelim ))
	#fout.write('\n')

	for i in range(len(methodList)) :
		fout.write('\n')
		for j in range(len(sampList)) :
			fout.write('{}{}'.format( resultsPR[i,j], textDelim ))
		fout.write('{}'.format(methodList[i]))
	#end loop
#end with



#TODO: get gene stats ... ??


print('\nDone.\n')