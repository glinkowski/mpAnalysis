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
dDir = 'pred04-test01'
dRoot = '../Dropbox/mp/output/'

# # Input names & locations
# useNtwk = 1		# network & samples to use (0 means fake)
# if useNtwk == 0 :
# 	eName = 'fakeNtwk01_g3e4t1'
# 	ePath = 'networks/'
# 	dRoot = 'outputFake/'
# else :
# 	eName = 'all_v3beta_g2e9t0'
# 	ePath = '../Dropbox/mp/networks/'
# 	dRoot = '../Dropbox/mp/output/'
# #end if

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
with open(dPath + 'parameters.txt', 'r') as fin :
	line = fin.readline()

	line = fin.readline()
	line = line.rstrip()
	lv = line.split(textDelim)
	eName = lv[1]

	line = fin.readline()
	line = line.rstrip()
	lv = line.split(textDelim)
	ePath = lv[1]
#end with

print("Creating the gene-index dictionary.")
geneDict = mp.readGenesFile(ePath, eName)
geneList = list(geneDict.keys())
geneList.sort()

print("Reading in the path names.")
pathDict = mp.readKeyFile(ePath, eName)
pathList = mp.removeInvertedPaths(pathDict)
del pathDict



# 2) Get the list of unique samples
dSubDirs = mp.getSubDirectoryList(dRoot+dDir)
sampSet = list()
prevSN = ''
numFolds = 0
prevNumFolds = 0
for si in dSubDirs :

	# extract sample name from folder name
	sv = si.split('/')
	sDir = sv[-2]
	sdv = sDir.split('-')
	sn = sdv[1]

	# count how many folds there are
	if sn != prevSN :
		if prevNumFolds != numFolds :
			print("ERROR: inconsistent number of folders per sample.")
			sys.exit
		#end if
		prevNumFolds = numFolds
		numFolds = 0
	else :
		numFolds += 1
#end loop
sampList = list(sampSet)
sampList.sort()



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
elif rp < 1 :
	print("ERROR: No ranked_paths... files exist in {}".format(dSubDirs[0]))
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
		methodList.append(fv[1])
#end loop
methodList.sort()

resultsROC = np.zeros( (len(methodList), len(dSubDirs)), dtype=np.float64)
resultsPR = np.zeros( (len(methodList), len(dSubDirs)), dtype=np.float64)
resultsAvgROC = np.zeros( (len(methodList), len(sampList)), dtype=np.float64)
resultsAvgPR = np.zeros( (len(methodList), len(sampList)), dtype=np.float64)
resultsMaxROC = np.zeros( (len(methodList), len(sampList)), dtype=np.float64)
resultsMaxPR = np.zeros( (len(methodList), len(sampList)), dtype=np.float64)
resultsPaths = np.zeros( (len(pathList), len(dSubDirs)), dtype=np.float64)
resultsGenes = np.zeros( (len(geneList), len(dSubDirs)), dtype=np.float64)
print("matrix sizes: {}, {}, {}, {}".format( len(resultsROC), len(resultsPR),
	len(resultsPaths), len(resultsGenes) ))



# 4) Create the Area Under the Curve tables

col = -1
#for si in dSubDirs[0:1] :
for si in dSubDirs :
	col += 1

	# Get data relating to each method
	row = -1
	for m in methodList :
		row += 1

		fn = 'ranked_genes-' + m + '.txt'
		FPR, recall, precision, numHid = vl.getAUCstats(si, fn)

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
		outName = si + 'AUC-' + m + '.png'
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
		sv = si.split('/')
		sdir = sv[-2]
		sdv = sdir.split('-')
		plt.suptitle(sdv[1]+', concealed = {}'.format(numHid)+
			', ROC area = {:.3}'.format(areaROC))
		plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.4, hspace=None)

		# Save the figure
		plt.savefig(outName)
		plt.close()
#end loop

# Get average AUC for each sample, across all folds
for i in range(len(sampList)) :
	left = i * numFolds
	right = left + numFolds

	resultsAvgROC[i] = np.mean(resultsROC[:,left:right], axis=0)
	resultsAvgPR[i] = np.mean(resultsPR[:,left:right], axis=0)
#end loop


#TODO: write this into a func

# Write the AUC tables to file(s)
with open(dDir + 'results-AUC_ROC.txt', 'w') as fout :
	fout.write('Area Under ROC Curve, per-sample')
	fout.write('\nnetwork:{}{}'.format(textDelim, eName))
	fout.write('\nfolds:{}{}'.format(textDelim, numFolds))
	fout.write('\n')

	for j in range(len(sampList)) :
		fout.write('{}{}'.format( sampList[j], textDelim ))
	#fout.write('\n')

	for i in range(len(methodList)) :
		fout.write('\n')
		for j in range(len(sampList)) :
			fout.write('{}{}'.format( resultsAvgROC[i,j], textDelim ))
		fout.write('{}'.format(methodList[i]))
	#end loop
#end with



#TODO: get path ranks per sample
#TODO: get gene stats ... ??