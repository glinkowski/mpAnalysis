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



####### ####### ####### ####### 
# PARAMETERS

# Variables/Quantities
percHide = [0, 10, 25, 33, 50]		# percent of genes to conceal


# Input names & locations
useNtwk = 0		# network & samples to use (0 means fake)

if useNtwk == 0 :
#	eName = 'fakeNtwk00_g2e3t10'
	eName = 'fakeNtwk01_g3e4t1'
	ePath = 'networks/'
	sPath = 'samplesFake/'
	oRoot = '../Dropbox/mp/outputFake/'
else :
	eName = 'all_v1_g2e11t0'
	ePath = '../Dropbox/mp/networks/'
	sPath = '../Dropbox/mp/samples-test1/'
	oRoot = '../Dropbox/mp/output/'
#end if

# Output path
oDirPrefix = 'pred03-batch'


# Output file extension (.txt vs .gz)
fileExt = '.gz'
#TODO: add if statement below to handle writing .gz vs .txt


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
fOutputName = 'parameters.gz'
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
		oSubDirList.append(oSubDir)

		# Read genes from sample & partition into test/train
		gAll = mp.readSampleFiles(sPath+s, True, True)
		gKnown, gHidden = mp.partitionSample(ePath, eName,
			oSubDir, gAll, p)

		if newVerbose :
#			print "Analyzing metapaths in sample: {}".format(s)
			print ( "  partitioned into {} known".format(len(gKnown)) +
				" and {} concealed genes ...".format(len(gHidden)) )
		#end if

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

#TODO: Break out into separate function
# need to pass: pathDict (can get pathList), oSubDirList, geneDict, oPath
#	? anything else ?
# TODO: pack get mxsize into this func

textDelim = '\t'
matrixDT = np.float32

Pxx = np.zeros([len(oSubDirList), len(pathList)], dtype=matrixDT)
Pyy = np.zeros([len(geneDict), len(pathList)], dtype=matrixDT)
Pxy = np.zeros([len(geneDict), len(pathList), len(oSubDirList)], dtype=matrixDT)
Sxy = np.zeros([len(geneDict), len(geneDict)], dtype=matrixDT)
SxySum = np.zeros([len(geneDict), len(pathList), len(oSubDirList)], dtype=matrixDT)
#pathScores = np.zeros([len(pathList), len(sNames)])
count = 0
for pi in xrange(len(pathList)) :
	count += 1

	# Load a metapath matrix into memory
	matrix = mp.getPathMatrix(pathDict[pathList[pi]], ePath, eName, mxSize)

	# Pyy is the diagonal of the path matrix
	Pyy[:,pi] = matrix.diagonal()


	# Sxy is the original PathSim metric b/t individual genes
	Sxy = matrix.transpose()
	Sxy = np.add( matrix, Sxy )
	PyyPxx = np.reshape( Pyy[:,pi], [1,Pyy[:,pi].shape[0]])
	PyyPxx = np.add( Pyy[:,pi], PyyPxx )
	PyyPxx = np.add( PyyPxx, 1 )
	Sxy = np.divide( Sxy, PyyPxx )

	# Calculate the path counts & sums
	for si in xrange(len(oSubDirList)) :
		sRows = matrix[oSampLists[si],:]
		sCols = matrix[:,oSampLists[si]]

		# Pxx is the number of connections of the sample to itself
		Pxx[si,pi] = sRows[:,oSampLists[si]].sum()

		# Pxy is the number of connections of sample to/from a given gene
		Pxy[:,pi,si] = np.add( sRows.sum(axis=0), sCols.sum(axis=1) )

		# SxySum is PathSim summed over the set X
		SxyCols = Sxy[:,oSampLists[si]]
		SxySum[:,pi,si] = SxyCols.sum(axis=1)

	#end loop

	if newVerbose :
		print("  Examined path {}, {}".format(count, pathList[pi]))
		print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))
	elif not (count % 25) :
		print("  Examined {} of {} paths...".format(count, len(pathList)))
		print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))
#end loop
print("Finished examining metapaths.")
print("    --elapsed time: {:.3} (s)".format(time.time()-tstart))


# Save the path counts to file(s)
tstart2 = time.time()

# Write Pyy to file
fname1 = 'Pyy'+fileExt
fout1 = gzip.open(oPath+fname1, 'wb')
#fout1 = open(oPath+fname1, 'wb')
fname2 = 'Pyy_mod'+fileExt
fout2 = gzip.open(oPath+fname2, 'wb')
#fout2 = open(oPath+fname2, 'wb')
#TODO: PyyMod = add(Pyy, 1)
#	Then discard Pyy, remove +1 in later steps

firstRow = True
for r in xrange(Pyy.shape[0]) :
	if firstRow :
		firstRow = False
	else :
		fout1.write('\n')
		fout2.write('\n')
	#end if

	firstCol = True
	for c in xrange(Pyy.shape[1]) :
		if firstCol :
			firstCol = False
		else :
			fout1.write('{}'.format(textDelim))
			fout2.write('{}'.format(textDelim))
		#end if

		fout1.write('{}'.format(Pyy[r,c]))
		fout2.write('{}'.format( (Pyy[r,c] + 1) ))
#end loop
fout1.close()
fout2.close()

# write Pxx, Pxy, SxySum; one for each sample (si)
for si in xrange(len(oSubDirList)) :

	path = oSubDirList[si]

	# Write Pxx to file
	fname = 'Pxx'+fileExt
	fout = gzip.open(oPath+fname, 'wb')

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
	fname1 = 'Pxy'+fileExt
	fout1 = gzip.open(path+fname1, 'wb')
#	fout1 = open(path+fname1, 'wb')
	fname2 = 'Pxy_mod'+fileExt
#	fout2 = open(path+fname2, 'wb')
	fout2 = gzip.open(path+fname2, 'wb')
	fname3 = 'Pxy_norm'+fileExt
	fout3 = gzip.open(path+fname3, 'wb')
#	fout3 = open(path+fname3, 'wb')
	fname4 = 'Pxy_mod_norm'+fileExt
	fout4 = gzip.open(path+fname4, 'wb')
#	fout4 = open(path+fname4, 'wb')
	fname5 = 'SxySum'+fileExt
	fout5 = gzip.open(path+fname5, 'wb')
#	fout5 = open(path+fname5, 'wb')

	PxyMod = np.divide( Pxy[:,:,si], np.add( Pyy, 1 ) )
	PxyColMax = np.add( np.amax(Pxy[:,:,si], axis=0), 0.0001 )
#	PxyColMax = np.multiply( PxyColMax, 100 )
	PxyColMax = np.divide( 100, PxyColMax )
	PxyModColMax = np.add( np.amax(PxyMod, axis=0), 0.0001 )
#	PxyModColMax = np.multiply( PxyModColMax, 100 )
	PxyModColMax = np.divide( 100, PxyModColMax )

	firstRow = True
	for r in xrange(Pxy.shape[0]) :
		if firstRow :
			firstRow = False
		else :
			fout1.write('\n')
			fout2.write('\n')
			fout3.write('\n')
			fout4.write('\n')
			fout5.write('\n')
		#end if

		firstCol = True
		for c in xrange(Pxy.shape[1]) :
			if firstCol :
				firstCol = False
			else :
				fout1.write('{}'.format(textDelim))
				fout2.write('{}'.format(textDelim))
				fout3.write('{}'.format(textDelim))
				fout4.write('{}'.format(textDelim))
				fout5.write('{}'.format(textDelim))
			#end if

			fout1.write('{}'.format(Pxy[r,c,si]))
			fout2.write('{}'.format( PxyMod[r,c] ))
#TODO: If I don't want to do this entry-wise, and I don't want to make
#	more matricies, I can at least to this row-wise:
#			fout3.write('{}'.format( int(round( Pxy[r,c,si] / PxyColMax[c] )) ))
#			fout4.write('{}'.format( int(round( PxyMod[r,c] / PxyModColMax[c] )) ))
			fout3.write('{}'.format( int(round( Pxy[r,c,si] * PxyColMax[c] )) ))
			fout4.write('{}'.format( int(round( PxyMod[r,c] * PxyModColMax[c] )) ))
			fout5.write('{}'.format( SxySum[r,c,si] ))
	#end loop
	fout1.close()
	fout2.close()
	fout3.close()
	fout4.close()
	fout5.close()
#end loop
print("Finished writing path counts.")
print("    --this step: {:.3} (s)".format(time.time()-tstart2))
print("    --total elapsed time: {:.3} (s)".format(time.time()-tstart))



print("\nDone.\n")