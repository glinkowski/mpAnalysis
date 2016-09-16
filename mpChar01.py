# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Characterization: define unique types of connections in set
#	Step 1:
#		create z-score features
#
# For each sample, perform once on the whole set and then
#	once on each of 4 folds. Partition set into Known Pos
#	and Hidden Pos.
# Then, create the metapath z-score feature, a per-gene
#	feature vector where each entry is the sum of that
#	gene's similarity to each member of the Known Pos set.
# ---------------------------------------------------------

import mpCharLib as cl
import time
import numpy as np
import gzip



######## ######## ######## ######## 
# PARAMETERS

# Input names & locations
useNtwk = 0		# network & samples to use (0 means fake)
if useNtwk == 0 :
#	eName = 'fakeNtwk00_g2e3t10'
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
	sPath = 'samplesFake/'
	oRoot = 'outputFake/'
else :
	eName = 'all_v3beta_g2e9t0'
	ePath = '../Dropbox/mp/networks/'
	sPath = '../Dropbox/mp/samples/MSIGDB/select01/'
#	sPath = '../Dropbox/mp/samples-test1/'
	oRoot = '../Dropbox/mp/output/'
#end if


# verbose feedback ?
verboseOutput = True

######## ######## ######## ######## 




######## ######## ######## ######## 
# PRIMARY FUNCTION

# Partition the samples and create z-score features
# Input ----
#	eName, str: folder containing network files
#	ePath, str: path to the network folder
#	sDir, str: directory containing samples
#	oDor, str: output directory
#	printFlag, bool: whether to enable terminal output
# Returns ----
#	nothing
# Creates ----
#	directory in oDir containing folder for each sample
#	(whole + 4 folds), containing gene partitions, and
#	a matrix of z-score feature vectors for that sample
def createFeatureZScore(eName, ePath, sDir, oRoot, printFlag) :

	######## ######## ######## ######## 
	# Parameters
	percHide = 0.25
#TODO: remove any superfluous params

	#1) define an appropriate output directory
	#	save parameters file
	#2) Read in the samples
	#	for each, create 4 folds
	#	Save samples & folds to oDir
	#3) Read in network info (gene names, etc)
	#4) for each mp & set, create z-score feature
	#5) save the feature matrix to each set



	# 1) Name & create a folder to store output files
	oSubDir = cl.nameOutputPath(oRoot, 'char01-batch')
	if printFlag :
		print("Files will be saved to {}".format(oSubDir))
	oDir = oRoot + oSubDir

	# Save experiment parameters to file
	fOutput = list()
	fOutput.append( ['date', 'network', 'ntwk path', '% hidden', 'samples'] )
	fOutput.append( [time.strftime('%d/%m/%Y'), eName, ePath, percHide, sDir] )
	fOutputName = 'parameters.txt'
	cl.writeGenericLists(oDir, fOutputName, fOutput)


	# 2a) Load the gene-index dict
	if printFlag :
		print("Creating the gene-index dictionary.")
	geneDict = cl.getGeneDictionary(ePath, eName)

	# 2b) Get the list of available paths
	if printFlag :
		print("Checking what paths are available ...")
	pathDict = cl.getPathDictionary(ePath, eName)
	pathList = cl.removeInvertedPaths(pathDict)


#TODO: this should be a function
	# Get expected matrix size
	print("Finding the matrix dimensions ...")
	fnMatrixZPad = 5
	#TODO: pack this into a function
	fname = (ePath+eName+"_MetaPaths/" +
		"{}.gz".format(str(0).zfill(fnMatrixZPad)) )
	mxSize = 0
	with gzip.open(fname, 'rb') as fin :
		for line in fin :
			mxSize += 1
	#end with



	# 3) Read & partition samples, save to output dir
	sNames = cl.getSampleNamesFromFolder(sDir)

	oSampLists = list()
	oSubDirList = list()

	# Read samples & create cross-validation folds
	sCount = 0
	for s in sNames :
		sCount += 1
		if printFlag :
			print("Collecting sample: {}, {}".format(sCount, s))

		# Read in genes from the full sample
		# Remove any genes not in the network
		# Convert sample names to indices
		gAll = cl.readSampleFiles(sDir + s, True, True)
		gAllValid, gIgnored = cl.checkListAgainstDictKeys(
			gAll, geneDict)
		giAll = [geneDict[g] for g in gAllValid]

		# Append the full sample to the lists
		oSampLists.append(giAll)
		oSubDir = '{}full-{}/'.format(oDir, s)
		oSubDirList.append( oSubDir )

		# Write the genes to file
		saveListToText(oSubDir, 'known.txt', gAllValid)
		saveListToText(oSubDir, 'hidden.txt', list() )
		saveListToText(oSubDir, 'ignored.txt', gIgnored)

		# Create the 4 folds
		random.shuffle(giAll)
		numFold = int(len(giAll) / 4)
		for i in range(4) :
			start = i * numFold
			stop = (i * numFold) + numFold
			gHidden = gAllValid[start:stop]
			gHidden.sort()
			gKnown = gAllValid[0:start]
			gKnown.extend( gAllValid[stop:0])
			gKnown.sort()

			# Append this fold to the lists
			giKnown = [geneDict[g] for g in gKnown]
			oSampLists.append(giKnown)
			oSubDirList.append( '{}part-{}-{02d}/'.format(sDir, s, i) )

			# Write the genes to file
			saveListToText(oSubDir, 'known.txt', gKnown)
			saveListToText(oSubDir, 'hidden.txt', gHidden)
			saveListToText(oSubDir, 'ignored.txt', list() )
		#end loop range(4)
	#end loop sNames



	# 4) Create the z-score features


	#	Build the feature vector matrices for each sample
	gFeatures = np.zeros( (len(geneDict), len(pathList),
		len(oSampLists)), dtype=np.float32 )

	# populate dimension 2 from each path
	dim2 = -1
	for p in pathList :
		dim2 += 1
		
		tpath = time.time()

		# load the path count matrix
		countMatrix = cl.getPathMatrix(pathDict[p], ePath, eName, mxSize)

		# Calculate z-score over the columns
		countAvg = np.mean(countMatrix, axis=0)
		countAvg = countAvg.reshape( (1, len(countAvg)) )
		countStD = np.std(countMatrix, axis=0)
		countStD = countStD.reshape( (1, len(countStD)) )
		countStD = np.add(countStD, 0.0001)

		# convert to z-scores on a row-by-row basis
		simMatrix = np.subtract(countMatrix, countAvg)
		simMatrix = np.divide(countMatrix, countStD)

		# zero out the diagonal (remove gene self-similarity)
		np.fill_diagonal(simMatrix, 0)

		# populate dimension 3 from each sample
		dim3 = -1
		for giList in oSampLists :
			dim3 += 1

			# calculate feature for this sample variant
			simSet = np.sum( simMatrix[:,giList], axis=1 )
			gFeatures[:,dim2,dim3] = simSet[:]
		#end loop

		if printFlag and not (dim2 % 25) :
			print("  Examined {} of {} paths...".format(dim2, len(pathList)))
			print("    --time per path: {:.3} (s)".format(time.time()-tpath))
#TODO: convert into a "time remaining / ETA"
	#end loop
	if printFlag :
		print("Finished examining matrix similarity matrices.")


#TODO: clean this
	# 5) ...
	Save the raw feature vectors to the sample folders
	i = -1
	twrite = time.time()
	for sDir in oSubDirList :
		i += 1

		mp.saveMatrixNumpy(gFeatures[:,:,i], 'features_ZScoreSim',
			sDir, False)
	#end loop
	print("Finished writing Z-Score Similarity feature vector files.")
	print("    --elapsed time: {:.3} (s)".format(time.time()-twrite))






#end def ######## ######## ######## 




######## ######## ######## ######## 
# FUNCTION CALL

tstart = time.time()

print("\nCreating z-score features for samples in {}".format(sPath))

# Main Function call
print("Calling function createFeatureZScore() ...")
#eDir = cl.concatenatePaths(ePath, eName)
createFeatureZScore(eName, ePath, sPath, oRoot, verboseOutput)

print("--elapsed time: {:.3} (s)".format(time.time()-tstart))



print("\nDone.\n")