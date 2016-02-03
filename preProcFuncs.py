# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# functions used in Pre-Processing of the network
#
# These functions were created to aid in the
#	pre-processing of a network, prior to calculated the
#	contained metapaths.
#
# Functions provided:
#	readKeepFile(fname)
#	readEdgeFile(datafile)
#	applyCorrections(edges, fname)
#	applyNormalization(edges, lowBound)
#	applyThreshold(edges, threshold)
#	applyKeepEdges(edges, kEdges)
#	applyKeepLists(edges, lGenes, kEdges, iEdges)
#	createNodeLists(edges, aGenes)
#	createModEdgeFileName(name, kEdges, kGenes, tHold)
#	writeModEdgeFilePlus(path, oname, nDict, gList, eArray)
#	createGeneMapping(gList)
#	createCountDict(aList)
#	createMatrixList(eArray, kEdges, iEdges, gList,	nDict)
#	saveMatrixText(matrix, mname, mpath, integer)
#	saveMatrixNumpy(matrix, mname, mpath)
#	clearFilesInDirectory(path)
#	saveMatrixList(mList, mNames, mGenes, mpath)
#	saveMatrixListPlus(mList, mNames, mGenes, mpath)
#	saveKeyFile(mDict, path)
#	saveGeneFile(mGenes, path)
#	createMPLengthOne(pList, pNames, path)
#	createMPLengthTwo(pList, pNames, path)
#	createMPLengthThree(pList, pNames, path)
#	createMPLengthFour(pList, pNames, path)
#	createMetaPaths(pList, pNames, gList, depth, path)
#	readPrimaryMatrices(nName, nPath)
# ---------------------------------------------------------

import os.path
#import os
import sys
import numpy as np
import re



####### ####### ####### ####### 
# PARAMETERS

# Data type used when loading edge file to memory:
nodeDT = np.dtype('a30')
# Whether to use the data-type for the matrices:
speedVsMemory = True	# True favors speed, disables dtype
# Data-type for the path matrices:
matrixDT = np.uint16
warnDTvalue = 65000
# Length to pad the matrix file names:
keyZPad = 5
# Whether to save uncompressed text version of matrix:
saveText = False
# Considering consecutive edges of same type
keepDouble = True
keepTriple = True
# File extension to use when saving the matrix
matrixExt = '.gz'	# '.txt' or '.gz' (gz is compressed)

####### ####### ####### ####### 



######## ######## ######## ########
# Function: Read in the keep file corresponding to
#	the specified network edge file
#	The keep file specifies which genes & edges to
#	keep in the final product, and provides the
#	necessary regex characters to identify them.
# Input:
#	fname, str - path & name to keep file
# Returns:
#	keepGenes - list of Regex expressions for the
#		genes to keep specified in file
#	keepEdges - list of Regex expressions for the
#		edges to keep specified in file
#	indirEdges - list of indirect edge types
def readKeepFile(fname) :

	# Make sure the keep file exists
	if not os.path.isfile(fname) :
		print "Please create a keep file: ", fname
		sys.exit()
	#end if

	# The lists to return
	keepGenes = list()
	loseGenes = list()
	keepEdges = list()
	indirEdges = list()

	# Read the file
	f = open(fname, "rb")
	line = f.readline()    # throw away the first line
#	for line in f :
	while line :
		# split the line by columns
		line = f.readline()
		line = line.rstrip()
		lv = line.split('\t')

		# Read in the gene types (# given in file)
		if lv[0] == 'GENE TYPES' :
			count = int(lv[2])
			for i in range(0, count) :
				line = f.readline()
				line = line.rstrip()
				lv = line.split('\t')
				if lv[2] == 'yes' :
					keepGenes.append(lv[1])
				else :
					loseGenes.append(lv[1])
				#end if
			#end loop
		# Read in the edge types (# given in file)
		elif lv[0] == 'EDGE TYPES' :
			count = int(lv[2])
			for i in range(0, count) :
				line = f.readline()
				line = line.rstrip()
				lv = line.split('\t')
				if lv[2] == 'yes' :
					keepEdges.append(lv[1])
				#end if
			#end loop
		elif lv[0] == 'INDIRECT' :
			count = int(lv[2])
			for i in range(0, count) :
				line = f.readline()
				line = line.rstrip()
				lv = line.split('\t')
				if lv[2] == 'yes' :
					indirEdges.append(lv[1])
			#end loop
		elif lv[0] == 'THRESHOLD' :
			tHold = float(lv[2])
		#end if
	#end loop

	return keepGenes, loseGenes, keepEdges, indirEdges, tHold

#end def ######## ######## ########



######## ######## ######## ########
# Function: Read in the Knowledge Graph
# Input:
#	fname, str - path & name to the network edge file 
# Returns:
#   Edges - (Nx4) matrix of char strings
#       each row is: node, node, edge weight, edge type
#   Nodes - dictionary of nodes in the edge list
#       key = node name
#       value = list of indices of the rows in
#           the edge matrix where node (key) appears
# matrix of Edges (N,4), set of Vertices
def readEdgeFile(datafile) :

	# get the number of lines in the file
	nLines = sum( 1 for line in open(datafile, "rb") )

	# assign space for edge list
#    dt = np.dtype('a30')
	dt = nodeDT
#    dt = np.dtype(str)
#    dt = np.dtype(object)

#    Edges = np.empty( [nLines,4], dtype=dt)
	Edges = np.empty( (nLines,4), dtype=dt)

	# dictionary to hold Node indices
	Nodes = dict()
	nodeSet = set()

	# Start reading from the file
	df = open(datafile, "rb")

	i = 0
	for line in df:
		# extract the data from the file
		line = line.rstrip()
		lv = line.split('\t')

		# insert into the edge list
		Edges[i,0] = lv[0]
		Edges[i,1] = lv[1]
		Edges[i,2] = lv[2]
		Edges[i,3] = lv[3]

		# add node locations to dict
		if (lv[0] in nodeSet) :
			Nodes[lv[0]].append(i)
		else :
			Nodes[lv[0]] = list()
			Nodes[lv[0]].append(i)
			nodeSet.add(lv[0])
		#end if
		if (lv[1] in nodeSet) :
			Nodes[lv[1]].append(i)
		else :
			Nodes[lv[1]] = list()
			Nodes[lv[1]].append(i)
			nodeSet.add(lv[1])
		#end if

		i += 1
	#end loop

	# close the data file
	df.close()

#    print "  file contained {:,} lines".format(nLines)
	return Edges, Nodes
#end def ######## ######## ########



######## ######## ######## ########
# Function: Fix known spelling mistakes, as outlined
#	in the corrections file. The file contains two
#	columns: the typo, the correction. If no file,
#	then leave the array unchanged.
# Input:
#	edges, array (N,4) - edge array
#		col 0: node, 1: node, 2: weight, 3: edge type 
# Returns:
#	nothing
#	Makes in-place corrections to the array
def applyCorrections(edges, fname) :

	# Make no changes if no file found
	if not os.path.isfile(fname) :
		print "No file found, no corrections made."
		return #edges
	#end if

	checkDict = dict()	# column 1
	checkSet = set()
	fixList = list()	# column 2

	cf = open(fname, "rb")
	index = 0
	for line in cf :
		line = line.rstrip()
		lv = line.split('\t')

		checkSet.add(lv[0])
		checkDict[lv[0]] = index
		index += 1
		fixList.append(lv[1])
	#end loop

#	print checkSet
#	print edges.shape

	for i in range(0, edges.shape[0]) :
		if edges[i,0] in checkSet :
			edges[i,0] = fixList[checkDict[edges[i,0]]]
		#end if
		if edges[i,1] in checkSet :
			edges[i,1] = fixList[checkDict[edges[i,1]]]
		#end if
		if edges[i,3] in checkSet :
#			print edges[i,3], fixList[checkDict[edges[i,3]]]
			edges[i,3] = fixList[checkDict[edges[i,3]]]
		#end if
	#end loop

	return #edges

#end def ######## ######## ########



######## ######## ######## ########
# Function: For each edge type, collect the weights
#	and normalize to [0,1].
# Input:
#	edges, array (N,4) - edge array
#		col 0: node, 1: node, 2: weight, 3: edge type 
#	lowBound, int - indicates what should be used as
#		the lower bound on the weight. If 0, then set
#		vmin = 0. Else, vmin = lowest weight.
# Returns:
#	nothing
#	Makes in-place corrections to the array
def applyNormalization(edges, lowBound) :

#TODO: check & fix lower bound condition. Don't want edges
#	to be given a weight of 0. There should be a small
#	weight applied even to the lowest-weighted edge

	# get the unique edge types in the network
	eTypes = list( np.unique(edges[:,3]) )
	eTypes.sort()
#	print "contains the following edges:"#.format(infile)
#	print eTypes


#	print "Normalizing edge type ..."
	for et in eTypes :
#		print "    {}".format(et)

		# extract the edge weights for only this edge type
		weightStr = edges[edges[:,3]==et, 2]
		# convert entries from str to float
		weightFlt = np.zeros((len(weightStr),1))
		for i in range(0,len(weightStr)) :
			weightFlt[i] = float(weightStr[i])
		#end loop

		if lowBound == 0 :
			vmin = 0
		else :
			vmin = float(np.amin(weightFlt))
		#end if

		vmax = float(np.amax(weightFlt))
		vdif = vmax - vmin
	#    print vmin, vmax, vdif

		for i in range(0, edges.shape[0]) :
			# calculate normalized value & update array
			if edges[i,3] == et :
				if (vdif == 0) :
					# if vmin=vmax, then there's only one value
					temp = 1
				else :
					# else change value to range [0,1]
					temp = (
						float(edges[i,2]) - vmin) / vdif
				#end if
#				print temp, vmin, vdif
				edges[i,2] = str(temp)
			#end if
		#end loop
		
	#end loop

	return
#end def ######## ######## ########



######## ######## ######## ########
# Function: Examine the edge weights. Throw out edges
#	that are below the threshold value. (Only keep
#	those that are at or above.)
# Input:
#	edges, array (N,4) - edge array
#		col 0: node, 1: node, 2: weight, 3: edge type 
#	threshold, float - the value against which to test
#		edge weights; throw out edges that are below
# Returns:
#	newEdges, str array - the modified edge array
def applyThreshold(edges, threshold) :
	# Not using del[], because each del is O(n)

#	if threshold > 1 :
#		print ("\nWARNING: threshold value {}".format(threshold)
#			+ " is outside the normalized range.")
#		print ("    If network has been normalized, all"
#			+ " edges will be kept.")
#	#end if

	# Get the weights in the network, count how
	#	many are at/above the threshold value
	count = 0
	for i in range(0, edges.shape[0]) :
		if float(edges[i,2]) >= threshold :
			count += 1
		#end if
	#end loop

	# If no edges are above threshold value ...
	if count == 0 :
		print ("\nWARNING: No edge weights above "
			+ " {}, returning unaltered array.\n".format(threshold))
		return edges
	#end if

	# temporary array to hold the new matrix
	newEdges = np.empty([count,4], dtype=nodeDT)
#	print newEdges.shape, count

	# Read the weight of each edge, and add to
	#	temp array if at or above threshold
	index = 0
	for i in range(0, edges.shape[0]) :
		weight = float(edges[i,2])
		if weight >= threshold :
#			print i, index
#			print edges[i,:]
			newEdges[index] = edges[i]
			index += 1
#			newEdges.append(edges[i])
		#end if
	#end loop

#	print newEdges

	# Replace current network with new one
#	edges = newEdges
#	return

#TODO: If desired threshold value is 0, don't throw
# out any edges. 
	return newEdges
#end def ######## ######## ########



######## ######## ######## ########
# Function: Examine the edge weights. Throw out edges
#	that are below the threshold value. (Only keep
#	those that are at or above.)
# Input:
#	edges, str array (N,4) - edge array
#		col 0: node, 1: node, 2: weight, 3: edge type 
#	kGenes, str list - regex of genes to keep
#		each entry is just the first four chars of
#		the gene name
#	kEdges, str list - edge types (full) to keep
#	threshold, float - the value against which to test
#		edge weights; throw out edges that are below
# Returns:
#	newEdges, str array - the modified edge array
# NOTE: This assumes the network has been altered
#	to just gene-gene edges !!!
def applyKeepEdges(edges, kEdges) :

	keepIndex = list()
	kEdgeSet = set(kEdges)

	for i in range(0, edges.shape[0]) :

		if edges[i,3] in kEdgeSet :
			keepIndex.append(i)
		#end if

#		# Throw out non-kept edges
#		if edges[i,3] not in kEdgeSet :
#			# Skip this edge
#			continue
#		#end if
#		keepIndex.append(i)

	#end loop

	newEdges = edges[keepIndex, :]
	return newEdges
#end def ######## ######## ########



######## ######## ######## ########
# Function: 
# Input:
#	edges, str array (N,4) - edge array
#		col 0: node, 1: node, 2: weight, 3: edge type 
#	kGenes, str list - regex of genes to keep
#		each entry is just the first four chars of
#		the gene name
#	kEdges, str list - edge types (full) to keep
#	threshold, float - the value against which to test
#		edge weights; throw out edges that are below
# Returns:
#	newEdges, str array - the modified edge array
def applyKeepLists(edges, lGenes, kEdges, iEdges) :

	keepIndex = list()
	kEdgeSet = set(kEdges)

#	print "----- applyKeepLists --------"
#	print kEdges
#	print kEdgeSet
#	for s in kEdgeSet : print s

	for i in range(0, edges.shape[0]) :

		# Throw out non-kept edges
		if edges[i,3] not in kEdgeSet :
			# Skip this edge
#			print "drop1", edges[i,:]
			continue
		#end if

		# list of matches to be found
		m0 = list()
		m1 = list()
		# Check nodes for matches (column 1 & 2)
		for gt in lGenes :
			m0.append( re.match(gt, edges[i,0]) )
			m1.append( re.match(gt, edges[i,1]) )
		#end loop
		# Throw out genes that match the non-keep list
		# Check for any match with the non-keep list
		if any(match != None for match in m1) :
#			print "drop2", edges[i,:]
			# Skip this edge
			continue
		#ASSUMPTION: for indirect edges, col 0 contains
		#	a non-gene node
		elif edges[i,3] not in iEdges :
			if any(match != None for match in m0) :
#				print "drop3", edges[i,:]
				# Skip this edge
				continue
		#end if

		# Finally, if no objections
		#	add this to the list to keep
		keepIndex.append(i)
	#end loop

	newEdges = edges[keepIndex, :]
	return newEdges
#end def ######## ######## ######## 



######## ######## ######## ########
# Function: Create a node dictionary from the current
#	edge list.
# Input:
#	edges, str array (N,4) - edge array
#		col 0: node, 1: node, 2: weight, 3: edge type 
# Returns:
#	nodeDict - node dictionary
#		key: str, name of node
#		value: list of int, list of indices where
#			these nodes occur in the edge list
def createNodeLists(edges, aGenes) :

	nodeDict = dict()
	nodeSet = set()
	geneSet = set()

	for i in range(0, edges.shape[0]) :

		# Add the first node to the dictionary,
		#	using a set for look-up speed
		if edges[i,0] in nodeSet :
			nodeDict[edges[i,0]].append(i)
		else :
			nodeDict[edges[i,0]] = list()
			nodeDict[edges[i,0]].append(i)
			nodeSet.add(edges[i,0])
		#end if

		# Add the second node to the dictionary,
		#	using a set for look-up speed
		if edges[i,1] in nodeSet :
			nodeDict[edges[i,1]].append(i)
		else :
			nodeDict[edges[i,1]] = list()
			nodeDict[edges[i,1]].append(i)
			nodeSet.add(edges[i,1])
		#end if


		# list of matches to be found
		m0 = list()
		m1 = list()
		# Check nodes for matches (column 1 & 2)
		for gt in aGenes :
			m0.append( re.match(gt, edges[i,0]) )
			m1.append( re.match(gt, edges[i,1]) )
		#end loop
		# Matches mean node is a gene; add to set
		if any(match != None for match in m0) :
			geneSet.add(edges[i,0])
		if any(match != None for match in m1) :
			geneSet.add(edges[i,1])
		#end if

	#end loop

	geneList = list(geneSet)
	geneList.sort()
	return nodeDict, geneList
#end def ######## ######## ######## 



######## ######## ######## ########
# Function: creates the name to use when saving the file
# Input:
#	name, str - the name of the original edge file
#	kEdges, str list - list of edge types kept
#	kGenes, str list - list of gene types kept
#	tHold, float - threshold value for kept edge weights
# Returns:
#	oname, str - the new network name for the new files
def createModEdgeFileName(name, kEdges, kGenes, tHold) :

	# save edge list, node dict, genes?
#	oname = name + "_g{1}t{0}_".format( int(tHold*100),
#		len(kGenes))

#	oname = ename + "_t{0}%_g{1}_".format( int(thresh*100),
#	    len(keepGenes))

#	for et in kEdges :
#		oname = oname + et[0]
#	#end loop

	oname = name + "_g{}e{}t{}".format( len(kGenes),
		len(kEdges), int(tHold*100) )
	return oname
#end def ######## ######## ######## 



######## ######## ######## ########
# Function: write the modified network to a file
#	This includes the nodeDict and geneList
# Input:
#	path, str - path to the folder to save the file
#	oname, str - new network name
#	nDict, dict - 
#		keys, node names as strings
#		values, int list of indexes: rows in eArray where
#			the node appears
#		eArray, (Nx4) array - the network
# Returns:
def writeModEdgeFilePlus(path, oname, nDict, gList, eArray) :

	newPath = path + oname + '/'

	# If folder doesn't exist, create it
	if not os.path.exists(newPath) :
		os.makedirs(newPath)
	#end if

	gfile = "genes.txt"
	nfile = "indices.txt"
	gf = open(newPath + gfile, "wb")
	nf = open(newPath + nfile, "wb")
	first = True
	for gene in gList :

		# Start every new line but the first with \n
		#   This avoids ending the file with a blank line
		if first == True :
			first = False
		else :
			gf.write("\n")
			nf.write("\n")
		#end if

		gf.write("{}".format(gene))

		nf.write("{}\t".format(gene, nDict[gene]))

		fIndex = True
		for item in nDict[gene] :
			if fIndex == True :
				fIndex = False
			else :
				nf.write(",")
			#end if
			nf.write("{}".format(item))
		#end loop

	#end loop

	gf.close()
	nf.close()

	ofile = "edge.txt"
	of = open(newPath + ofile, "wb")

	first = True
	for i in range(0, eArray.shape[0]) :
		if first == True :
			first = False
		else :
			of.write("\n")
		#end if

		of.write("{}\t{}\t{}\t{}".format(eArray[i,0],
			eArray[i,1], eArray[i,2], eArray[i,3]))

	#end loop

	of.close()

	return
#end def ######## ######## ######## 



######## ######## ######## ########
# Function: create a mapping of genes to matrix indices
#	ASSUMPTION: gList is properly ordered to match matrix
# Input:
#	gList, str list: ordered list of gene names
# Returns:
#	gDict, dict
#		key, str: gene names
#		value, int: row/col where gene appears in matrix
def createGeneMapping(gList) :

	gDict = dict()
	numG = 0
	for gene in gList :
		gDict[gene] = numG
		numG += 1
	#end loop

	return gDict
#end def ######## ######## ######## 



######## ######## ######## ########
# Function: create a dict of node counts
# Input:
# Returns:
def createCountDict(aList) :
	# From the list of genes in this network, dictionary
	#	will return a count of how many times each appears

	aSet = set()
	aDict = dict()

	for item in aList :

		if item in aSet :
			aDict[item] += 1
		else :
			aSet.add(item)
			aDict[item] = 1
		#end if
	#end loop

	return aDict
#end def ######## ######## ######## 
######## ######## ######## ########
# Function: create a list of the primary matrices
# Input:
#	path, str - path to the folder to save the file
#	oname, str - new network name
#	nDict, dict - 
#		keys, node names as strings
#		values, int list of indexes: rows in eArray where
#			the node appears
#		eArray, (Nx4) array - the network
# Returns:
def createMatrixList(eArray, kEdges, iEdges, gList,
	nDict): #, path, oname) :

	# data-type for the path matrices
	dt = matrixDT
	warnVal = warnDTvalue

	# How far from Std Dev will be kept as "Medium"
	stdGap = 0.75

	# Get a gene-to-index mapping to use with the
	#	newly-created matrices
	gDict = createGeneMapping(gList)
	numG = len(gDict.keys())

	# Create a set for faster membership checks
	gSet = set(gList)

	mList = list()
	mNames = list()

	iEdges.sort()
	kEdges.sort()

	# Define cut-offs for indirect edges
	# ASSUMPTION: indirect edges, col 0 is non-gene node
	iNames = dict()
		# key: edge type + sm/md/lg
		# value: integer tuple w/ low & high cutoffs

	for et in iEdges :

		# Get the degree of each node for this edge type
		nodeList = list(eArray[eArray[:,3]==et, 0])
		nodeCount = createCountDict(nodeList)

		# Get the mean & std dev for the node degrees
		mean = int(round(np.mean(nodeCount.values())))
		std = int(round(np.std(nodeCount.values())))
		stdPart = int(round(std * stdGap))

		# Create the breakdown by small/med/large
		tSml = [max(0, mean-(2*std)), mean - stdPart]
		tMed = [tSml[1]+1, mean + stdPart]
		tLrg = [tMed[1]+1, mean+(2*std)]

		# Save those tuples to the look-up dict
		iNames[et+"_sm"] = tSml
		iNames[et+"_md"] = tMed
		iNames[et+"_lg"] = tLrg
	#end loop

#	print iNames

	# Initialize empty file
	# This file stores what types were kept & how many edges
#	fn = path + oname + 'types.txt'
#	fet = open(fn, 'wb')
#	#fet.write("{}{}{}\n".format(et, delim, count))
#	fet.close()

	# Start creating matrices
	for et in kEdges :

		# Remove indirect edges if needed
		if et in iEdges :

			# Create lists corresponding to the new edge types
			term_sm = list()
			term_md = list()
			term_lg = list()

			# Separate out the edges of type et
			thisArray = eArray[eArray[:,3]==et]

			# Get the count of each node with that edge
			checkSet = set()
			termDict = dict()
			for row in thisArray :

				# Only add if the node is not a gene
				if row[0] not in gSet :
					# Else, increment count by 1
					if row[0] in checkSet :
						termDict[row[0]] += 1
#						print "if", row[0]
					# If not yet added, then set = 1
					else :
#						print "else", row[0]
						checkSet.add(row[0])
						termDict[row[0]] = 1
					#end if
				#end if

				## Only add if the node is not a gene
				#if row[1] not in gSet :
				#	# If not yet added, then set = 1
				#	if row[1] in checkSet :
				#		termDict[row[1]] == 1
				#		checkSet.add(row[1])
				#	# Else, increment count by 1
				#	else :
				#		termDict[row[1]] += 1
				##end if
			#end loop

			# Assign to groups according to node degree
			for term in termDict.keys() :

				# get the tuples to check
				tSml = iNames[row[3] + "_sm"]
				tMed = iNames[row[3] + "_md"]
				tLrg = iNames[row[3] + "_lg"]


				if (tSml[0] <= termDict[term] <= tSml[1]) :
					term_sm.append(term)
				elif (tMed[0] <= termDict[term] <= tMed[1]) :
					term_md.append(term)
				elif (tLrg[0] <= termDict[term] <= tLrg[1]) : 
					term_lg.append(term)
				#end if
			#end loop


			# Create the first (small) matrix
#			print "Creating matrix from edge type: {}".format(et+' < '+str(cutf[et][1]))
			if speedVsMemory :
				thisM = np.zeros([numG, numG])
			else :
				thisM = np.zeros([numG,numG], dtype=dt)
			#end if
			count = 0
			for term in term_sm :
				# list of all edges with this term
				rowList = nDict[term]
				# Join all genes connected to term X
				for i in range(0, len(rowList)) :
					for j in range(0+1, len(rowList)) :
						# Find two genes joined by term
						# ASSUME: terms always in col 0, genes in col 1
						gA = eArray[i,1]
						gB = eArray[j,1]
						# Increment the entry(s) in the array (symmetric)
						thisM[gDict[gA],gDict[gB]] += 1
						thisM[gDict[gB],gDict[gA]] += 1
						count += 1
					#end loop
				#end loop
			#end loop			
			if (count > 0) :
				# save to a file
#				fn = mpath + oname + et + '_sm'
#				print "    saving to {}".format(fn)
#				np.save(fn, thisM)
#				# This file stores what types were kept & how many edges
#				fn = mpath + oname + 'types.txt'
#				fet = open(fn, 'ab')
#				fet.write("{}\t{}\n".format(et+'_sm', count))
#				fet.close()

				# ERROR CHECK: verify counts fit within
				#	specified data type
				if np.amax(thisM) > warnVal :
					print ("WARNING: Path counts exceed" +
						"{}, change data-type.".format(warnVal))
				#end if

				mList.append(thisM)
				mNames.append(et+"_sm")
			#end if

			# Create the second (medium) matrix
#			print "Creating matrix from edge type: {}".format(et+' < '+str(cutf[et][2]))
			if speedVsMemory :
				thisM = np.zeros([numG, numG])
			else :
				thisM = np.zeros([numG,numG], dtype=dt)
			#end if
			count = 0
			for term in term_md :
				# list of all edges with this term
				rowList = nDict[term]
				# Join all genes connected to term X
				for i in range(0, len(rowList)) :
					for j in range(0+1, len(rowList)) :
						# Find two genes joined by term
						# ASSUME: terms always in col 0, genes in col 1
						gA = eArray[i,1]
						gB = eArray[j,1]
						# Increment the entry(s) in the array (symmetric)
						thisM[gDict[gA],gDict[gB]] += 1
						thisM[gDict[gB],gDict[gA]] += 1
						count += 1
					#end loop
				#end loop
			#end loop
			if (count > 0) :
				# save to a file
#				fn = mpath + oname + et + '_md'
#				print "    saving to {}".format(fn)
#				np.save(fn, thisM)
#				# This file stores what types were kept & how many edges
#				fn = mpath + oname + 'types.txt'
#				fet = open(fn, 'ab')
#				fet.write("{}\t{}\n".format(et+'_md', count))
#				fet.close()

				# ERROR CHECK: verify counts fit within
				#	specified data type
				if np.amax(thisM) > warnVal :
					print ("WARNING: Path counts exceed" +
						"{}, change data-type.".format(warnVal))
				#end if

				mList.append(thisM)
				mNames.append(et+"_md")
			#end if
			
			# Create the third (large) matrix
#			print "Creating matrix from edge type: {}".format(et+' < '+str(cutf[et][3]))
			if speedVsMemory :
				thisM = np.zeros([numG, numG])
			else :
				thisM = np.zeros([numG,numG], dtype=dt)
			#end if
			count = 0
			for term in term_lg :
				# list of all edges with this term
				rowList = nDict[term]
				# Join all genes connected to term X
				for i in range(0, len(rowList)) :
					for j in range(0+1, len(rowList)) :
						# Find two genes joined by term
						# ASSUME: terms always in col 0, genes in col 1
						gA = eArray[i,1]
						gB = eArray[j,1]
						# Increment the entry(s) in the array (symmetric)
						thisM[gDict[gA],gDict[gB]] += 1
						thisM[gDict[gB],gDict[gA]] += 1
						count += 1
					#end loop
				#end loop
			#end loop
			if (count > 0) :
				# save to a file
#				fn = mpath + oname + et + '_lg'
#				print "    saving to {}".format(fn)
#				np.save(fn, thisM)
#				# This file stores what types were kept & how many edges
#				fn = mpath + oname + 'types.txt'
#				fet = open(fn, 'ab')
#				fet.write("{}\t{}\n".format(et+'_lg', count))
#				fet.close()

				# ERROR CHECK: verify counts fit within
				#	specified data type
				if np.amax(thisM) > warnVal :
					print ("WARNING: Path counts exceed" +
						"{}, change data-type.".format(warnVal))
				#end if
				
				mList.append(thisM)
				mNames.append(et+"_lg")
			#end if
			



		# If already direct, create the matrix
		else :
			if speedVsMemory :
				thisM = np.zeros([numG, numG])
			else :
				thisM = np.zeros([numG,numG], dtype=dt)
			#end if
			count = 0

	#		print "Creating matrix from edge type: {}".format(et)
			thisArray = eArray[eArray[:,3]==et]
			# increment entry at (i,j) = (gene0,gene1)
			for row in thisArray :
				thisM[gDict[row[0]],gDict[row[1]]] += 1
				thisM[gDict[row[1]],gDict[row[0]]] += 1
				count += 1
			#end loop

			# This file stores what types were kept & how many edges
#			fn = mpath + oname + '.types.txt'
#			fet = open(fn, 'ab')
#			fet.write("{}\t{}\n".format(et, count))
#			fet.close()

			# ERROR CHECK: verify counts fit within
			#	specified data type
			if np.amax(thisM) > warnVal :
				print ("WARNING: Path counts exceed" +
					"{}, change data-type.".format(warnVal))
			#end if
				
			mList.append(thisM)
			mNames.append(et)
		#end if
	#end loop

	return mList, mNames
#end def ######## ######## ######## 



######## ######## ######## ########
# Function: save the given matrix as a .txt file
# Input:
#	matrix, (NxN) list - the values to save
#	mname, str - name of the file to save
#	mpath, str - path to the folder to save the file
#	integer, bool - True means save values as int()
# Returns:
#	nothing
def saveMatrixText(matrix, mname, mpath, integer) :

	# If folder doesn't exist, create it
	if not os.path.exists(mpath) :
		os.makedirs(mpath)
	#end if

	# Open the file
	fout = open(mpath + mname + ".txt", "wb")

	# Write to the file
	firstR = True
	for i in range(0, matrix.shape[0]) :

		# if not the first row, start with \n
		if firstR :
			firstR = False
		else :
			fout.write("\n")
		#end if

		firstC = True
		for j in range(0, matrix.shape[1]) :

			# if not the first col, start with \t
			if firstC :
				firstC = False
			else :
				fout.write("\t")
			#end if

			# Write the value to file
			#	If integer = True, write as an integer
			if integer :
				fout.write("{}".format( int(matrix[i,j]) ))
			else :
				fout.write("{}".format( matrix[i,j] ))
			#end if

		#end loop
	#end loop

	fout.close()

	return
#end def ######## ######## ########



######## ######## ######## ########
# Function: save the given matrix as a .npy file
# Input:
#	matrix, (NxN) list - the values to save
#	mname, str - name of the file to save
#	mpath, str - path to the folder to save the file
#	integer, bool - True means save values as int()
# Returns:
#	nothing
def saveMatrixNumpy(matrix, mname, mpath) :

	# If folder doesn't exist, create it
	if not os.path.exists(mpath) :
		os.makedirs(mpath)
	#end if

	# Write to the file
#	np.save(mpath+mname, matrix)
	np.savetxt(mpath + mname + matrixExt, matrix, fmt='%u')
#NOTE: In this case, the text file from savetxt() is much
#	smaller than the binary file from save()

	#ERROR CHECK: also save a text version of the matrix
	if saveText :
		saveMatrixText(matrix, mname, mpath, True)
	#end if

	return
#end def ######## ######## ######## 



######## ######## ######## ########
# Function: delete the files within a directory
# Input:
#	path, str - the directory to delete
# Returns:
#	nothing
# NOTE: This only works if no sub-folders
def clearFilesInDirectory(path) :

	# Get the list of files in the directory
	filelist = ([f for f in os.listdir(path)
		if os.path.isfile(os.path.join(path, f))])

	# Remove all the files in that folder
	for f in filelist :
		os.remove(os.path.join(path, f))

	# Delete the folder
#	os.rmdir(path)

	return
#end def ######## ######## ######## 



######## ######## ######## ########
# Function: save a list of matrices
# Input:
#	mList, list of NxN matrices - the matrices to save
#	mNames, list of str - names of the paths in matrix
#	mGenes, list of str - names of genes in the matrix
#	mpath, str - path to the folder to save the file
# Returns:
#	nothing
# Creates:
#	
def saveMatrixList(mList, mNames, mGenes, mpath) :

	# Number of digits to zero-pad the file name/number
#	zpad = 5
	zpad = keyZPad

	# If folder doesn't exist, create it
	if not os.path.exists(mpath) :
		os.makedirs(mpath)
	# otherwise, delete files in the directory
	else :
		clearFilesInDirectory(mpath)
	#end if

	# This file gives the corresponding gene names for
	#	each row/col of the matrix (rows & cols are same)
	fgene = open(mpath+"genes.txt", "wb")
	firstline = True
	for gene in mGenes :
		if firstline :
			firstline = False
		else :
			fgene.write("\n")
		#end if
		fgene.write("{}".format(gene))
	#end if

	# This file tells which matrix corresponds to which path
	fkey = open(mpath+"key.txt", "wb")

	num = 0
	firstline = True
	for i in range(0, len(mNames)) :

		# Write to the legend file
		if firstline :
			firstline = False
		else :
			fkey.write("\n")
		#end if
		fkey.write("{:05d}\t{}".format(num, mNames[i]))

		# Save each matrix as the corresponding number
		saveMatrixNumpy(mList[i], str(num).zfill(zpad),
			mpath)

		# VERIFICATION: save as a text-readable file
		saveMatrixText(mList[i], "t"+str(num).zfill(zpad),
			mpath, True)

		num += 1
	#end loop

	return
#end def ######## ######## ######## 



######## ######## ######## ########
# Function: save a list of matrices
# Input:
#	mList, list of NxN matrices - the matrices to save
#	mNames, dict
#		key, str: metapath names
#		value, int: corresponding index number for mList
#	mGenes, list of str - names of genes in the matrix
#	mpath, str - path to the folder to save the file
# Returns:
#	nothing
# Creates:
#	
def saveMatrixListPlus(mList, mNames, mGenes, mpath) :

	# Number of digits to zero-pad the file name/number
#	zpad = 5
	zpad = keyZPad

	# If folder doesn't exist, create it
	if not os.path.exists(mpath) :
		os.makedirs(mpath)
	#end if

	# This file gives the corresponding gene names for
	#	each row/col of the matrix (rows & cols are same)
	fgene = open(mpath+"genes.txt", "wb")
	firstline = True
	for gene in mGenes :
		if firstline :
			firstline = False
		else :
			fgene.write("\n")
		#end if
		fgene.write("{}".format(gene))
	#end if
	fgene.close()


	# Get the sorted list of all paths
	nameList = list(mNames.keys())
	nameList.sort()

	# This file tells which matrix corresponds to which path
	fkey = open(mpath+"key.txt", "wb")
	fkey.write("NOTE: 't' - use matrix transpose\n")
	firstline = True
	for name in nameList :
		if firstline :
			firstline = False
		else :
			fkey.write("\n")
		#end if
#		fkey.write("{:05d}\t{}".format(mNames[name], name))
#		fkey.write("{:05d}".format(mNames[name][0]))
		fkey.write("{}".format( str(mNames[name][0]).zfill(zpad) ))
		if mNames[name][1] == True :
			fkey.write(",t")
		else :
			fkey.write(", ")
		fkey.write("\t{}".format(name))
	#end loop
	fkey.close()

	for i in range(0, len(mList)) :
		# Save each matrix as the corresponding number
		saveMatrixNumpy(mList[i], str(i).zfill(zpad),
			mpath)

		# VERIFICATION: save as a text-readable file
		saveMatrixText(mList[i], "t"+str(i).zfill(zpad),
			mpath, True)
	#end loop

	return
#end def ######## ######## ######## 



######## ######## ######## ########
# Function: save the key file for the metapath matrices
# Input:
#	mDict, dict
#		key, str: metapath names
#		value, [int, : corresponding index number for mList
#				bool] : True means use matrix transpose
#	path, str - path to the folder to save the file
# Returns: nothing
# Creates: a legend, mapping the path type on the right
#	to the path matrix file on the left, where 't'
#	indicates the transpose of that matrix should be used
def saveKeyFile(mDict, path) :

	# The length to pad the file name / matrix number
	zpad = keyZPad

	# If folder doesn't exist, create it
	if not os.path.exists(path) :
		os.makedirs(path)
	#end if

	# Get the sorted list of all paths
	nameList = list(mDict.keys())
	nameList.sort()

	# This file tells which matrix corresponds to which path
	fkey = open(path+"key.txt", "wb")
	fkey.write("NOTE: 't' means use matrix transpose\n")
	firstline = True
	for name in nameList :
		if firstline :
			firstline = False
		else :
			fkey.write("\n")
		#end if
#		fkey.write("{:05d}\t{}".format(mNames[name], name))
#		fkey.write("{:05d}".format(mNames[name][0]))
		fkey.write("{}".format( str(mDict[name][0]).zfill(zpad) ))
		if mDict[name][1] == True :
			fkey.write(",t")
		else :
			fkey.write(", ")
		fkey.write("\t{}".format(name))
	#end loop
	fkey.close()

	return
#end def ######## ######## ######## 
######## ######## ######## ########
# Function: save the list of genes
# Input:
#	mGenes, list of str - list of genes in matrix
#		ASSUMPTION: list is already properly ordered
#	path, str - path to the folder to save the file
# Returns: nothing
# Creates: file containing ordered list of genes to
#	use as the row/col headers for path matrices
def saveGeneFile(mGenes, path) :

	# If folder doesn't exist, create it
	if not os.path.exists(path) :
		os.makedirs(path)
	#end if

	# This file gives the corresponding gene names for
	#	each row/col of the matrix (rows & cols are same)
	fgene = open(path+"genes.txt", "wb")
	firstline = True
	for gene in mGenes :
		if firstline :
			firstline = False
		else :
			fgene.write("\n")
		#end if
		fgene.write("{}".format(gene))
	#end if
	fgene.close()

	return
#end def ######## ######## ######## 



######## ######## ######## ######## 
# Function: Read in the key.txt file regarding the 
#	metapath matrices
# Input ----
#	path, str: path to the network files
#	name, str: name of the network to use
# Returns ----
#	keyDict, dict
#		key, str: name of metapath
#		value, tuple: int is matrix/file ID number
#			bool where True means use matrix transpose
def readKeyFilePP(path) :

	fname = path + "key.txt"
	# ERROR CHECK: verify file exists
	if not os.path.isfile(fname) :
		print ( "ERROR: Specified file doesn't exist:" +
			" {}".format(fname) )
		sys.exit()
	#end if

	# The item to return
	keyDict = dict()

	# Read in the file
	fk = open(fname, "rb")
	firstline = True
	for line in fk :
		# skip the first line
		if firstline :
			firstline = False
			continue
		#end if

		# separate the values
		line = line.rstrip()
		lk = line.split('\t')
		lv = lk[0].split(',')

		transpose = False
		if lv[1] == "t" :
			transpose = True
		#end if

		# add to the dict
		keyDict[lk[1]] = [int(lv[0]), transpose]
	#end loop
	fk.close()

	return keyDict
#end def ######## ######## ######## 


def createMPLengthOne(pList, pNames, path) :
	mNum = 0
	zpad = keyZPad
	mDict = dict()
	for i in range(0, len(pNames)) :
		saveMatrixNumpy(pList[i], str(mNum).zfill(zpad), path)
		mDict[pNames[i]] = [mNum, False]
		mNum += 1
	#end loop
	saveKeyFile(mDict, path)
	return
#end def ######## ######## ######## 
def createMPLengthTwo(pList, pNames, path) :
	mDict = readKeyFilePP(path)
	mNum = len(mDict)
	zpad = keyZPad

	for i in range(0, len(pNames)) :
		for j in range(i, len(pNames)) :
			# Optionally skipping consecutive edges
			if not keepDouble :
				if i==j :
					continue
			#end if

			# The name of this path
			name = pNames[i] + "-" + pNames[j]
			# The name of the reversed path
			nameRev = pNames[j] + "-" + pNames[i]

			# Create new matrix if file doesn't already exist
			if not os.path.isfile(path + str(mNum).zfill(zpad) + matrixExt) :
				newM = np.dot(pList[i], pList[j])
				saveMatrixNumpy(newM, str(mNum).zfill(zpad), path)
			#end if

			# Add the matrix name & number to mDict
			if i == j :
				# If name == nameRev (ie: typeA-typeA)
				# Then add this matrix to the list
				mDict[name] = [mNum, False]
			else :
				# Add this path & note the reverse path
				mDict[name] = [mNum, False]
				#	Reverse path uses transpose
				mDict[nameRev] = [mNum, True]
			#end if
			mNum += 1
		#end loop
	#end loop

	saveKeyFile(mDict, path)
	return
#end def ######## ######## ######## 
def createMPLengthThree(pList, pNames, path) :
	mDict = readKeyFilePP(path)
	mNum = len(mDict)+1
	zpad = keyZPad

	checkSet = set()
	for i in range(0, len(pNames)) :
		for j in range(0, len(pNames)) :
			# Optionally skipping consecutive edges
			if not keepDouble :
				if i==j :
					continue
			#end if
			for k in range(0, len(pNames)) :
				# Optionally skipping consecutive edges
				if not keepDouble :
					if j==k :
						continue
				#end if
				# Skip if i=j=k (three in a row)
				if not keepTriple :
					if (i==j) and (j==k) :
						continue
				#end if

#				print "        creating {}, {}-{}-{}".format((mNum+1), i,j,k)

				# The name of this path
				name = ( pNames[i] + "-" +
					pNames[j] + "-" + pNames[k] )
				# The name of the reversed path
				nameRev = ( pNames[k] + "-" +
					pNames[j] + "-" + pNames[i] )

				# Verify this path wasn't yet calculated
				#	if it has been, skip it
				if name not in checkSet :
					checkSet.add(name)

					# Create new matrix if file doesn't already exist
					if not os.path.isfile(path + str(mNum).zfill(zpad) + matrixExt) :
						# Calculate the matrix
						temp = np.dot(pList[i], pList[j])
						newM = np.dot(temp, pList[k])
						# Save the data
						saveMatrixNumpy(newM, str(mNum).zfill(zpad), path)
					#end if

					mDict[name] = [mNum, False]

					# Check the reverse path (the transpose)
					if nameRev not in checkSet :
						checkSet.add(nameRev)
						# Save the data
						saveMatrixNumpy(newM, str(mNum).zfill(zpad), path)
						mDict[nameRev] = [mNum, True]
					#end if
				#end if

				mNum += 1
			#end loop
		#end loop
	#end loop

	saveKeyFile(mDict, path)
	return
#end def ######## ######## ######## 
def createMPLengthFour(pList, pNames, path) :
	mDict = readKeyFilePP(path)
	mNum = len(mDict)
	zpad = keyZPad

	checkSet = set()
	for h in range(0, len(pNames)) :
		for i in range(0, len(pNames)) :
			# Optionally skipping consecutive edges
			if not keepDouble :
				if h==i :
					continue
			#end if
			for j in range(0, len(pNames)) :
				# Optionally skipping consecutive edges
				if not keepDouble :
					if i==j :
						continue
				#end if
				# Skip if h=i=j (three in a row)
				if not keepTriple :
					if (h==i) and (i==j) :
						continue
				#end if
				for k in range(0, len(pNames)) :
					# Optionally skipping consecutive edges
					if not keepDouble :
						if j==k :
							continue
					#end if
					# Skip if i=j=k (three in a row)
					if not keepTriple :
						if (i==j) and (j==k) :
							continue
					#end if

					# The name of this path
					name = ( pNames[h] + "-" + pNames[i] +
						"-" + pNames[j] + "-" + pNames[k] )
					# The name of the reversed path
					nameRev = ( pNames[k] + "-" + pNames[j] +
						"-" + pNames[i] + "-" + pNames[h] )

					# Verify this path wasn't yet calculated
					#	if it has been, skip it
					if name not in checkSet :
						checkSet.add(name)

						# Create new matrix if file doesn't already exist
						if not os.path.isfile(path + str(mNum).zfill(zpad) + matrixExt) :
							# Calculate the matrix
							temp1 = np.dot(pList[h], pList[i])
							temp2 = np.dot(temp1, pList[j])
							newM = np.dot(temp2, pList[k])
							# Save the data
							saveMatrixNumpy(newM, str(mNum).zfill(zpad), path)
						#end if

						mDict[name] = [mNum, False]

						# Check the reverse path (the transpose)
						if nameRev not in checkSet :
							checkSet.add(nameRev)
							# Save the data
							saveMatrixNumpy(newM, str(mNum).zfill(zpad), path)
							mDict[nameRev] = [mNum, True]
						#end if
					#end if

					mNum += 1
				#end loop
			#end loop
		#end loop
	#end loop

	saveKeyFile(mDict, path)
	return
#end def ######## ######## ######## 
######## ######## ######## ########
# Function: save a list of matrices
# Input:
#	pList, list of NxN matrices - the primary matrices
#		ie: the 1-level path matrices
#	pNames, dict
#		key, str: metapath names
#		value, int: corresponding index number for mList
#	mGenes, list of str - names of genes in the matrix
#	mpath, str - path to the folder to save the file
# Returns: nothing
# Creates: The set of path matrices in the appropriate
#	folder. These are simply named numerically, with a
#	key/legend file provided. The list of genes used as
#	row/col headers is also saved to that folder.	
def createMetaPaths(pList, pNames, gList, depth, path) :

	maxDepth = 4
	if depth > maxDepth :
		print ( "WARNING: Can only calculate up to " +
			"{}-step metapaths.".format(maxDepth) )
	elif depth < 1 :
		print ( "WARNING: Requested metapaths of length" +
			" {};".format(depth) +
			" Will return only 1-step paths.")
		depth = 1
	#end if

	# Check if folder at specified path exists
	# If directory exists, emtpy it
	if os.path.exists(path) :
		clearFilesInDirectory(path)
	# Else, create the folder
	else :
		os.makedirs(path)
	#end if

	# The items to return
	mDict = dict()

	# The number of matrices created
	mNum = 0
	# The length to pad the file name / matrix number
	zpad = keyZPad

	#ERROR CHECK: verify gene list & matrix dimensions
	if len(gList) != pList[0].shape[0] :
		print ( "ERROR: The provided list of genes" +
			" does not match the matrix. No paths created.")
		return
	elif pList[0].shape[0] != pList[0].shape[1] :
		print ( "ERROR: The primary path matrix passed" +
			" is not square.")
		return
	#end if

	# Save the list of genes to file
	saveGeneFile(gList, path)

	#-------------------------------
	# Create the 1-step paths
	createMPLengthOne(pList, pNames, path)

	if depth < 2 :
#		saveKeyFile(mDict, path)
#		return mList, mDict
		return
	#end if
	print "    finished creating paths of length 1"

	#-------------------------------
	## Create the 2-step paths
	createMPLengthTwo(pList, pNames, path)
	print "    finished creating paths of length 2"
	
	if depth < 3 :
#		saveKeyFile(mDict, path)
#		return mList, mDict
		return
	#end if

	#-------------------------------
	# Create the 3-step paths
	createMPLengthThree(pList, pNames, path)
	print "    finished creating paths of length 3"

	if depth < 4 :
#		saveKeyFile(mDict, path)
#		return mList, mDict
		return
	#end if

	#-------------------------------
	# Create the 4-step paths
	createMPLengthFour(pList, pNames, path)
	print "    finished creating paths of length 4"

#	return mList, mDict
	return
#end def ######## ######## ######## 



######## ######## ######## ########
# Function: read in the primary matrices
# Input: 
#	nName, str - name of the network
#	nPath, str - path to the network
# Returns: 
#	pList, list of NxN matrices - the primary matrices
#		ie: the 1-level path matrices
#	pNames, dict
#		key, str: metapath names
#		value, int: corresponding index number for mList
def readPrimaryMatrices(nName, nPath) :

	# Check for folder existence
	path = nPath + nName + "_Primaries/"
	if not os.path.exists(path) :
		print "ERROR: Path doesn't exist: {}".format(path)
		sys.exit()
	#end if

	# Items to return
	pNames = list()
	pList = list()

	# Read in the key file
	fn = open(path + "key.txt", "rb")
	for line in fn :
		line = line.rstrip()
		lv = line.split('\t')

		pNames.append(lv[1])

		if speedVsMemory :
			if os.path.isfile( path + lv[0] + '.gz' ) :
				pList.append( np.loadtxt(path + lv[0] + '.gz') )
			elif os.path.isfile( path + lv[0] + '.txt' ) :
				pList.append( np.loadtxt(path + lv[0] + '.txt') )
			else :
				print ("ERROR: Unknown file name and extension" +
					" for matrix {}.".format(lv[0]))
				sys.exit()
			#end if
		else :
			if os.path.isfile( path + lv[0] + '.gz' ) :
				pList.append( np.loadtxt(path + lv[0] + '.gz', dtype=matrixDT) )
			elif os.path.isfile( path + lv[0] + '.txt' ) :
				pList.append( np.loadtxt(path + lv[0] + '.txt', dtype=matrixDT) )
			else :
				print ("ERROR: Unknown file name and extension" +
					" for matrix {}.".format(lv[0]))
				sys.exit()
			#end if
		#end if
	#end loop

	return pNames, pList
#end def ######## ######## ######## 