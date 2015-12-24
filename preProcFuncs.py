

import os.path
import sys
import numpy as np
import re


nodeDT = np.dtype('a30')



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
				keepEdges.append(lv[1])
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
def applyKeepLists(edges, lGenes, kEdges) :

	keepIndex = list()
	kEdgeSet = set(kEdges)

#	print kEdges
#	print kEdgeSet
#	for s in kEdgeSet : print s

	for i in range(0, edges.shape[0]) :

		# Throw out non-kept edges
		if edges[i,3] not in kEdgeSet :
			# Skip this edge
#			print "drop", edges[i,3]
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
		if any(match != None for match in m0) :
			# Skip this edge
			continue
		elif any(match != None for match in m1) :
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





