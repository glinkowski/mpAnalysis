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
dDir = 'pred04-set04'
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

# Make new pathDict to give index of pathList
pathDict = dict()
idx = -1
for item in pathList :
	idx += 1
	pathDict[item] = idx
#end loop



# 2) Get the list of unique samples
dSubDirs = mp.getSubDirectoryList(dRoot+dDir)
sampSet = set()
prevSN = ''
numFolds = 1
prevNumFolds = 1
for si in dSubDirs :

	# extract sample name from folder name
	sv = si.split('/')
	sDir = sv[-2]
	sdv = sDir.split('-')
	sn = sdv[1]
	sampSet.add(sn)

	# count how many folds there are
	if sn != prevSN :
#		if prevNumFolds != numFolds :
#			print("ERROR: inconsistent number of folders per sample.")
#			sys.exit()
#		#end if
		prevNumFolds = numFolds
		numFolds = 1
	else :
		numFolds += 1
	#end if
	prevSN = sn
#end loop
if prevNumFolds != numFolds :
	print("ERROR: inconsistent number of folders per sample.")
	sys.exit()
#end if
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
#print("matrix sizes: {}, {}, {}, {}".format( resultsROC.shape, resultsAvgROC.shape, resultsPaths.shape, resultsGenes.shape ))



# 4) Create the Area Under the Curve tables
print("Finding AUCs for each sample ...")

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
		plt.suptitle(sdv[1]+'\n{}, concealed = {}'.format(m, numHid)+
			', ROC area = {:.3}'.format(areaROC))
		plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.4, hspace=None)

		# Save the figure
		plt.savefig(outName)
		plt.close()
#end loop


#print(numFolds)
#print("{}, {}".format( resultsROC.shape, resultsAvgROC.shape ))
#print (resultsROC)
# Get average AUC for each sample, across all folds
for i in range(len(sampList)) :
	left = i * numFolds
	right = left + numFolds

#	newCol = np.mean(resultsROC[:,left:right], axis=1)
#	print("L:{}, R:{}".format(left, right))
#	print(newCol)
#	resultsAvgROC[:,i] = newCol[:]
	resultsAvgROC[:,i] = np.mean(resultsROC[:,left:right], axis=1)
	resultsMaxROC[:,i] = np.amax(resultsROC[:,left:right], axis=1)
#	resultsAvgROC[i] = np.mean(resultsROC[:,left:right], axis=0)
	resultsAvgPR[:,i] = np.mean(resultsPR[:,left:right], axis=1)
	resultsMaxPR[:,i] = np.amax(resultsPR[:,left:right], axis=1)
#end loop



# 5) Output to file(s)

#TODO: write this into a func
#print (resultsAvgROC)

# # Write the AUC tables to file(s)
# with open(dPath + 'results-AUC_ROC.txt', 'w') as fout :
# 	fout.write('Area Under ROC Curve, per-sample')
# 	fout.write('\nnetwork:{}{}'.format(textDelim, eName))
# 	fout.write('\nfolds:{}{}'.format(textDelim, numFolds))
# 	fout.write('\n')

# 	for j in range(len(sampList)) :
# 		fout.write('{}{}'.format( sampList[j], textDelim ))
# 	#fout.write('\n')

# 	for i in range(len(methodList)) :
# 		fout.write('\n')
# 		for j in range(len(sampList)) :
# 			fout.write('{}{}'.format( resultsAvgROC[i,j], textDelim ))
# 		fout.write('{}'.format(methodList[i]))
# 	#end loop
# #end with

print("Writing the AUC tables (ROC & PR) ...")

# Write the AUC tables to file(s) for the per-sample results
fileNames = (['results-AUC_ROC_Avg.txt', 'results-AUC_ROC_Max.txt',
	'results-AUC_PR_Avg.txt', 'results-AUC_PR_Max.txt'])
fileHeaders = (['Average -- Area Under ROC Curve, per-sample', 'Maximum -- Area Under ROC Curve, per-sample',
	'Average -- Area Under PR Curve, per-sample', 'Maximum -- Area Under PR Curve, per-sample'])
fileData = ([resultsAvgROC, resultsMaxROC, resultsAvgPR, resultsMaxPR])
for f in range(4) :
	with open(dPath + fileNames[f], 'w') as fout :
		fout.write(fileHeaders[f])
		fout.write('\nnetwork:{}{}'.format(textDelim, eName))
		fout.write('\nfolds:{}{}'.format(textDelim, numFolds))
		fout.write('\n\n')

		for j in range(len(sampList)) :
			fout.write('{}{}'.format( sampList[j], textDelim ))
		#fout.write('\n')

		for i in range(len(methodList)) :
			fout.write('\n')
			for j in range(len(sampList)) :
				fout.write('{}{}'.format( fileData[f][i,j], textDelim ))
			fout.write('{}'.format(methodList[i]))
		#end loop
#end loop

# Write the AUC tables to file(s) for the COMPLETE results (every fold)
fileNames = (['results-AUC_ROC_All.txt', 'results-AUC_PR_All.txt'])
fileHeaders = (['Area Under ROC Curve, all folds', 'Area Under PR Curve, all folds'])
fileData = ([resultsROC, resultsPR])
for f in range(2) :
	with open(dPath + fileNames[f], 'w') as fout :
		fout.write(fileHeaders[f])
		fout.write('\nnetwork:{}{}'.format(textDelim, eName))
		fout.write('\nfolds:{}{}'.format(textDelim, numFolds))
		fout.write('\n\n')

		for j in range(len(dSubDirs)) :
			sv = dSubDirs[j].split('/')
			fout.write('{}{}'.format( sv[-2], textDelim ))
		#fout.write('\n')

		for i in range(len(methodList)) :
			fout.write('\n')
			for j in range(len(dSubDirs)) :
				fout.write('{}{}'.format( fileData[f][i,j], textDelim ))
			fout.write('{}'.format(methodList[i]))
		#end loop
#end loop



# 6) Create the path score tables
print("Finding path scores for each sample ...")

pScoreFolds = np.zeros( (len(pathList), len(dSubDirs)) )
pScoreSamples = np.zeros( (len(pathList), len(sampList)) )

fcol = -1
#for si in dSubDirs[0:1] :
for si in dSubDirs :
	fcol += 1

	pScoreEvery = np.zeros( (len(pathList), len(methodList)) )

	# Get data relating to each method
	ecol = -1
	for m in methodList :
		ecol += 1

		fn = 'ranked_paths-' + m + '.txt'
		# Skip the ranked_paths file if missing
		if not os.path.isfile(si + fn) :
			continue
		with open(si + fn, 'r') as fin :
			firstLine = fin.readline()
			for line in fin :
				line = line.rstrip()
				lv = line.split(textDelim)

				idx = pathDict[lv[1]]
				pScoreEvery[idx,ecol] = float(lv[0])
		#end with

		# normalize to [-1, 1]
		absVals = np.absolute( pScoreEvery[:,ecol] )
		divVal = np.amax( absVals )
		divVal = divVal + 0.0001	# hack so as not to divide by 0
		pScoreEvery[:,ecol] = np.divide(pScoreEvery[:,ecol], divVal)
#		print(np.amax(pScoreEvery[:,ecol]))

		# place average/sum (?) into larger table
		pScoreFolds[:,fcol] = np.sum(pScoreEvery, axis=1)
#end loop


# Get average path scores for each sample, across all folds
for i in range(len(sampList)) :
	left = i * numFolds
	right = left + numFolds

	pScoreSamples[:,i] = np.mean(pScoreFolds[:,left:right], axis=1)
#end loop



# 7) Output to file(s)


print("Writing the path score files ...")

with open(dPath + 'results-PathScores_Avg.txt', 'w') as fout :
	fout.write('Path Scores for each sample, averaged across folds')
	fout.write('\nnetwork:{}{}'.format(textDelim, eName))
	fout.write('\nfolds:{}{}'.format(textDelim, numFolds))
	fout.write('\n\n')

	for j in range(len(sampList)) :
		fout.write('{}{}'.format( sampList[j], textDelim ))
	#fout.write('\n')

	for i in range(len(pathList)) :
		fout.write('\n')
		for j in range(len(sampList)) :
			fout.write('{}{}'.format( pScoreSamples[i,j], textDelim ))
		fout.write('{}'.format(pathList[i]))
#end with


with open(dPath + 'results-PathScores_All.txt', 'w') as fout :
	fout.write('Path Scores for each fold, averaged across all methods')
	fout.write('\nnetwork:{}{}'.format(textDelim, eName))
	fout.write('\nfolds:{}{}'.format(textDelim, numFolds))
	fout.write('\n\n')

	for j in range(len(dSubDirs)) :
		sv = dSubDirs[j].split('/')
		fout.write('{}{}'.format( sv[-2], textDelim ))
	#fout.write('\n')

	for i in range(len(pathList)) :
		fout.write('\n')
		for j in range(len(dSubDirs)) :
			fout.write('{}{}'.format( pScoreFolds[i,j], textDelim ))
		fout.write('{}'.format(pathList[i]))
#end with





#TODO: get gene stats ... ??


print('\nDone.\n')