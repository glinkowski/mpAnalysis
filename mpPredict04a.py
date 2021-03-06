# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Approach 2: learn top paths from PathSim sum, find genes
#	Version 2, Step 1 
#
# The first approach tries to first find the most
#	important metapaths, then rank genes by similarity
#	along those. This approach instead finds the set
#	similarity along each metapath, then tries to learn
#	the most important paths.
#
# Step 1: Build a per-gene feature vector where each entry is
#	the sum of that gene's similarity to each member of
#	the set in question. In this case, the PathSim matrices
#	have already been created as part of the network
#	pre-processing.
# ---------------------------------------------------------

import mpLibrary as mp
import time
import numpy as np
import gzip



####### ####### ####### ####### 
# PARAMETERS

# Variables/Quantities
percHide = [25, 25, 25, 25, 25]  # 5 x 25%
	# percent of genes to conceal


# Input names & locations
useNtwk = 1		# network & samples to use (0 means fake)
if useNtwk == 0 :
#	eName = 'fakeNtwk00_g2e3t10'
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
	sPath = 'samplesFake/'
	oRoot = 'outputFake/'
else :
	eName = 'all_v3beta_g2e9t0'
	ePath = '../Dropbox/mp/networks/'
	sPath = '../Dropbox/mp/samplesAchilles1st/subset03/'
#	sPath = '../Dropbox/mp/samples-test1/'
	oRoot = '../Dropbox/mp/output/'
#end if

# Output path
oDirPrefix = 'pred04-batch'


# verbose feedback ?
newVerbose = True

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()
print("")

mp.setParamVerbose(newVerbose)



# 0) Name & create a folder to store output files
oDirectory = mp.nameOutputPath(oRoot, oDirPrefix)
print("Files will be saved to {}".format(oDirectory))
oPath = oRoot + oDirectory

# Save experiment parameters to file
fOutput = list()
fOutput.append( ['date', 'network', 'ntwk path', '% hidden', 'samples'] )
fOutput.append( [time.strftime('%d/%m/%Y'), eName, ePath, percHide, sPath] )
fOutputName = 'parameters.txt'
mp.writeGenericLists(oPath, fOutputName, fOutput)



# 1) Load the gene-index dict
print("Creating the gene-index dictionary.")
geneDict = mp.readGenesFile(ePath, eName)



# 2) Get list of all samples in folder
sNames = mp.getSampleNamesFromFolder(sPath)



# 3) Read the samples & create cross-validation folds
oSampLists = list()
oSubDirList = list()

if not newVerbose :
	print("Collecting & partitioning samples ...")
	print("  there are {} samples".format(len(sNames)))

index = 0
sCount = 1
for s in sNames :
	if newVerbose :
		print("Collecting sample: {}, {}".format(sCount, s))
	sCount += 1


	count = -1
	for p in percHide :

		# Need a folder for each sample set
		count += 1
		oSubDir = '{}{:03d}-{}-{:02d}/'.format(oPath, p,
			s, count)
		oSubDirList.append(oSubDir)

		# Read genes from sample & partition into test/train
		gAll = mp.readSampleFiles(sPath+s, True, True)
		gKnown, gHidden = mp.partitionSample(ePath, eName,
			oSubDir, gAll, p)

		if newVerbose :
			print ( "  partitioned into {} known".format(len(gKnown)) +
				" and {} concealed genes ...".format(len(gHidden)) )
		#end if

		# Convert sample into list of indices
		giKnown = mp.convertToIndices(gKnown, geneDict)
		oSampLists.append(giKnown)
#end loop
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))



# 4) Get the list of available paths
print("Checking what paths are available ...")
pathDict = mp.readKeyFile(ePath, eName)
pathList = mp.removeInvertedPaths(pathDict)
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))

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
print("  mxSize = {}".format(mxSize))
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))



# 5) Get the similarity measure from PathSim matrices
#	Build the feature vector matrices for each sample
gFeatures = np.zeros( (len(geneDict), len(pathList),
	len(oSampLists)), dtype=np.float32 )

# populate dimension 2 from each path
dim2 = -1
for p in pathList :
#for p in [pathList[219], pathList[220]] :
	dim2 += 1

#	print("{}: {}".format(pathDict[p], p))

	# load the PathSim matrix
	simMatrix = mp.getSimMatrix(pathDict[p], ePath,
		eName, mxSize)
	# zero out the diagonal
	# NOTE: want to measure each gene's similarity to the set,
	#	excluding itself
	np.fill_diagonal(simMatrix, 0)

	# populate dimension 3 from each sample
	dim3 = -1
	for giList in oSampLists :
		dim3 += 1

#TODO: Explore this thought...
#	Why are length-2 paths so dominant?
#	Should I remove the diagonal from the PathSim sum?
#	Does the diagonal over-weight b/c of similarity to self?
#		# get the diagonal, the similarity of each gene to self
#		simDiag = np.diag(simMatrix)
#
		# calculate feature for this sample variant
		simSet = np.sum( simMatrix[:,giList], axis=1 )
#
#		# remove the similarity of the gene to itself,
#		#	sum of similiarity over the set, minus similarity to self
#
		gFeatures[:,dim2,dim3] = simSet[:]
	#end loop

	if newVerbose :
		print("  Examined path {}, {}".format(dim2, p))
		print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))
	elif not (dim2 % 25) :
		print("  Examined {} of {} paths...".format(dim2, len(pathList)))
		print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))
#end loop
print("Finished examining matrix similarity matrices.")
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))



# 6) Save the raw feature vectors to the sample folders
i = -1
for sDir in oSubDirList :
	i += 1

	mp.saveMatrixNumpy(gFeatures[:,:,i], 'features_PathSim',
		sDir, False)
#end loop
print("Finished writing PathSim feature vector files.")
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))



print("\nDone.\n")