# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
#
# A different approach to ranking meta-paths.
# In this case, the collected counts between genes and the
#	training set will be normalized and used to weight
#	the metapaths by predictability, or importance. This
#	can be throught of as using feedback from the PathSim
#	metric to rank genes in order to influence and automate
#	the decision of which paths to use for prediction.
#
# This is the first part of the process: collecting data
#	from all the metapaths for a batch of samples. This
#	will also create multiple versions of each sample using
#	different percentages of hidden genes.
# ---------------------------------------------------------

import mpLibrary as mp
import time
import numpy as np
import gzip


# For testing & verification
#import random
#random.seed(42)



####### ####### ####### ####### 
# PARAMETERS

# Variables/Quantities
percHide = [0, 10, 25, 33, 50]		# percent of genes to conceal
#nRandSamp = 200		# number of random samples to compare

# options for Node Binning
#useBinning = True
#nodeBins = [0.3, .66]
#binType = 'all'


# Input names & locations
useNtwk = 1		# network & samples to use (0 means fake)
if useNtwk == 0 :
#	eName = 'fakeNtwk00_g2e3t10'
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
	sPath = 'samplesFake/'
else :
	eName = 'all_v1_g2e11t0'
	ePath = '../Dropbox/mp/networks/'
	sPath = '../Dropbox/mp/samples-test2/'
#end if

# Output path
#oRoot = 'outputFake/'
oRoot = '../Dropbox/mp/output/'
oDirPrefix = 'pred03-batch'


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
fOutput.append( ['date', 'network', '% hidden'] )
fOutput.append( [time.strftime('%d/%m/%Y'), eName, percHide] )
fOutputName = 'parameters.txt'
mp.writeGenericLists(oPath, fOutputName, fOutput)


# 1) Load the gene-index dict
print("Creating the gene-index dictionary.")
geneDict = mp.readGenesFile(ePath, eName)


# 2) Get list of all samples in folder
sNames = mp.getSampleNamesFromFolder(sPath)
print sNames


# 3) Read the samples & create cross-validation folds
oSampLists = list()
oSubDirList = list()

index = 0
for s in sNames :
	print("Analyzing metapaths in sample: {}".format(s))
	for p in percHide :

		# Need a folder for each sample set
		oSubDir = oPath+'{0:03d}'.format(p)+'-'+s+'/'
#		oSubDir = oPath+p.zfill(3)+s+'/'
		oSubDirList.append(oSubDir)

		# Read genes from sample & partition into test/train
		gAll = mp.readSampleFiles(sPath+s, True, True)
		gKnown, gHidden = mp.partitionSample(ePath, eName,
			oSubDir, gAll, p)

#		print "Analyzing metapaths in sample: {}".format(s)
#		print ( "  partitioned into {} known".format(len(gKnown)) +
#			" and {} concealed genes ...".format(len(gHidden)) )

		# Convert sample into list of indices
		gIndices = mp.convertToIndices(gKnown, geneDict)
		oSampLists.append(gIndices)

#end loop
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))


# 4) Get the list of available paths
print("Checking what paths are available ...")
pathDict = mp.readKeyFile(ePath, eName)
pathList = mp.removeInvertedPaths(pathDict)
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))

# Get expected matrix size
print "Finding the matrix dimensions ..."
fnMatrixZPad = 5
#TODO: pack this into a function
fname = (ePath+eName+"_MetaPaths/" +
	"{}.gz".format(str(0).zfill(fnMatrixZPad)) )
mxSize = 0
with gzip.open(fname, 'rb') as fin :
	for line in fin :
		mxSize += 1
#end with
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))



# 5) Calculate path counts 

textDelim = '\t'

Pxx = np.zeros([len(oSubDirList), len(pathList)])
Pyy = np.zeros([len(geneDict), len(pathList)])
Pxy = np.zeros([len(geneDict), len(pathList), len(oSubDirList)])
#pathScores = np.zeros([len(pathList), len(sNames)])
count = 0
for pi in xrange(len(pathList)) :
	count += 1

	# Load a metapath matrix into memory
	matrix = mp.getPathMatrix(pathDict[pathList[pi]], ePath, eName, mxSize)

	# Pyy is the diagonal of the path matrix
	Pyy[:,pi] = matrix.diagonal()
#	print "Pyy: ", Pyy

	# Calculate the path counts & sums
	for si in xrange(len(oSubDirList)) :
		sRows = matrix[oSampLists[si], :]
		sCols = matrix[:, oSampLists[si]]

		# Pxx is the number of connections of the sample to itself
		Pxx[si, pi] = sRows[:, oSampLists[si]].sum()
#		print "Pxx: ", Pxx

		# Pxy is the number of connections of sample to/from a given gene
		Pxy[:,pi,si] = sRows.sum(axis=0) + sCols.sum(axis=1)

	#end loop

	if newVerbose :
		print("  Examined path {}".format(pathList[pi]))
		print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))
	elif not (count % 25) :
		print("  Examined {} of {} paths...".format(count, len(pathList)))
		print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))
#end loop
print("Finished examining metapaths.")
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))


# Save the path counts to file(s)

# Write Pyy to file
fname = 'Pyy.txt'
fname2 = 'Pyy-mod.txt'
fout = open(oPath+fname, 'wb')
fout2 = open(oPath+fname2, 'wb')
firstRow = True
for r in xrange(Pyy.shape[0]) :
	if firstRow :
		firstRow = False
	else :
		fout.write('\n')
		fout2.write('\n')
	#end if

	firstCol = True
	for c in xrange(Pyy.shape[1]) :
		if firstCol :
			firstCol = False
		else :
			fout.write('{}'.format(textDelim))
			fout2.write('{}'.format(textDelim))
		#end if

		fout.write('{}'.format(Pyy[r,c]))
		fout2.write('{}'.format( (Pyy[r,c] + 1) ))
#end loop
fout.close()
fout2.close()


for si in xrange(len(oSubDirList)) :

	path = oSubDirList[si]

	# Write Pxx to file
	fname = 'Pxx.txt'
	fout = open(oPath+fname, 'wb')
	firstCol = True
	for c in xrange(Pxx.shape[1]) :
		if firstCol :
			firstCol = False
		else :
			fout.write('{}'.format(textDelim))
		#end if
		fout.write('{}'.format(Pxx[si,c]))
	#end loop
	fout.close()


	# Write Pxy to file, along with normed version(s)

	fname = 'Pxy.txt'
	fout = open(path+fname, 'wb')

	fname2 = 'Pxy-mod.txt'
	fout2 = open(path+fname2, 'wb')

#	PxyMod = np.divide(Pxy[:,:,si], np.add(Pyy, 1) )
	PxyMod = np.divide(Pxy[:,:,si], (Pyy + 1) )
#	print Pxy[0,0,si], PxyMod[0,0], (Pyy[0,0] + 1)
	PxyColMax = np.add( np.amax(Pxy[:,:,si], axis=0), 0.0001)
	PxyModColMax = np.add( np.amax(PxyMod, axis=0), 0.0001)
#	print Pxy.shape, PxyColMax.shape
#TODO: further vectorization

	fname3 = 'Pxy-norm.txt'
	fout3 = open(path+fname3, 'wb')

	fname4 = 'Pxy-mod-norm.txt'
	fout4 = open(path+fname4, 'wb')

	firstRow = True
	for r in xrange(Pxy.shape[0]) :
		if firstRow :
			firstRow = False
		else :
			fout.write('\n')
			fout2.write('\n')
			fout3.write('\n')
			fout4.write('\n')
		#end if

		firstCol = True
		for c in xrange(Pxy.shape[1]) :
			if firstCol :
				firstCol = False
			else :
				fout.write('{}'.format(textDelim))
				fout2.write('{}'.format(textDelim))
				fout3.write('{}'.format(textDelim))
				fout4.write('{}'.format(textDelim))
			#end if

			fout.write('{}'.format(Pxy[r,c,si]))
#			fout2.write('{}'.format( (Pxy[r,c,si] / (Pyy[r,c] + 1)) ))
			fout2.write('{}'.format( PxyMod[r,c] ))
			fout3.write('{}'.format( int(round(Pxy[r,c,si] / PxyColMax[c] * 100 ) )))
			fout4.write('{}'.format( int(round(PxyMod[r,c] / PxyModColMax[c] * 100 ) )))
	#end loop
	fout.close()
	fout2.close()
	fout3.close()
	fout4.close()
#end loop
print("Finished writing path counts.")
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))




print("\nDone.\n")