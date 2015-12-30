

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
	oname = name + "_g{1}t{0}_".format( int(tHold*100),
		len(kGenes))
#	oname = ename + "_t{0}%_g{1}_".format( int(thresh*100),
#	    len(keepGenes))
	for et in kEdges :
		oname = oname + et[0]
	#end loop

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

	gfile = oname + ".genes.txt"
	nfile = oname + ".indices.txt"
	gf = open(path + gfile, "wb")
	nf = open(path + nfile, "wb")
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

	ofile = oname + ".edges.txt"
	of = open(path + ofile, "wb")

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
# Input:
# Returns:
def createGeneMapping(gList) :
	# From the list of genes in this network, dictionary
	#	will indicate what index to use in matrix
	# ASSUMPTION: file is sorted

	gDict = dict()
	numG = 0
	for gene in gList :
		gDict[gene] = numG
		numG += 1
	#end loop

	return gDict

#end def ######## ######## ######## ######## ######## ######## ########
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
def createMatrixList(eArray, kEdges, iEdges, gList, nDict, path, oname) :

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

	# where to save the matrix files
	mpath = path + "matrices/"

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

#	iNames = dict()
#	for et in iEdges :
#		valList = list()
#
##		eValStr = eArray[eArray[:,3] == et][:,2]
#		eValStr = eArray[eArray[:,3] == et, 2]
#
#		print "-----inside createMatrixList----------"
##		print eArray
##		print eArray[eArray[:,3] == et, :]
##		print et, eValStr
#
#		# Convert the strings to float
#		eVals = list()
#		eVals = [float(x) for x in eValStr]
##		print eVals
#
#		mean = np.mean(eVals)
#		std = np.std(eVals)
#
#		tSml = [max(0, mean-(2*std)), mean-(std*stdGap)]
#		tMed = [tSml[1]+1, mean + (std*stdGap)]
#		tLrg = [tMed[1]+1, mean+(2*std)]
#
##		iNames.append(e+"_sm", tSml)
##		iNames.append(e+"_md", tMed)
##		iNames.append(e+"_lg", tLrg)
#
#		iNames[et+"_sm"] = tSml
#		iNames[et+"_md"] = tMed
#		iNames[et+"_lg"] = tLrg
#	#end loop
#
#	print iNames


	# Initialize empty file
	# This file stores what types were kept & how many edges
	fn = path + oname + 'types.txt'
	fet = open(fn, 'wb')
	#fet.write("{}{}{}\n".format(et, delim, count))
	fet.close()

	# Start creating matrices
	for et in kEdges :

#		print et, kEdges, iEdges

		# Remove indirect edges if needed
		if et in iEdges :

#			print et, iEdges

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

#				print row
#				print termDict

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

#			print term_sm, term_md, term_lg


			# Create the first (small) matrix
#			print "Creating matrix from edge type: {}".format(et+' < '+str(cutf[et][1]))
			thisM = np.zeros([numG,numG])
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
				fn = mpath + oname + et + '_sm'
				print "    saving to {}".format(fn)
				np.save(fn, thisM)
				# This file stores what types were kept & how many edges
				fn = mpath + oname + 'types.txt'
				fet = open(fn, 'ab')
				fet.write("{}\t{}\n".format(et+'_sm', count))
				fet.close()

				mList.append(thisM)
				mNames.append(et+"_sm")
			#end if

			# Create the second (medium) matrix
#			print "Creating matrix from edge type: {}".format(et+' < '+str(cutf[et][2]))
			thisM = np.zeros([numG,numG])
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
				fn = mpath + oname + et + '_md'
				print "    saving to {}".format(fn)
				np.save(fn, thisM)
				# This file stores what types were kept & how many edges
				fn = mpath + oname + 'types.txt'
				fet = open(fn, 'ab')
				fet.write("{}\t{}\n".format(et+'_md', count))
				fet.close()

				mList.append(thisM)
				mNames.append(et+"_md")
			#end if
			
			# Create the third (large) matrix
#			print "Creating matrix from edge type: {}".format(et+' < '+str(cutf[et][3]))
			thisM = np.zeros([numG,numG])
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
				fn = mpath + oname + et + '_lg'
				print "    saving to {}".format(fn)
				np.save(fn, thisM)
				# This file stores what types were kept & how many edges
				fn = mpath + oname + 'types.txt'
				fet = open(fn, 'ab')
				fet.write("{}\t{}\n".format(et+'_lg', count))
				fet.close()

				mList.append(thisM)
				mNames.append(et+"_lg")
			#end if
			



		# If already direct, create the matrix
		else :
			thisM = np.zeros([numG,numG])
			count = 0

	#		print "Creating matrix from edge type: {}".format(et)
			thisArray = eArray[eArray[:,3]==et]
			# increment entry at (i,j) = (gene0,gene1)
			for row in thisArray :
				thisM[gDict[row[0]],gDict[row[1]]] += 1
				thisM[gDict[row[1]],gDict[row[0]]] += 1
				count += 1
			#end loop

			# save to a file
			fn = mpath + oname + "." + et
			print "    saving to {}".format(fn)
			np.save(fn, thisM)

			#ERROR CHECK: save to a text file
			fn = mpath + oname + "." + et + '.txt'
			print "    saving to {}".format(fn)
			np.savetxt(fn, thisM, delimiter='\t')

			# This file stores what types were kept & how many edges
			fn = mpath + oname + '.types.txt'
			fet = open(fn, 'ab')
			fet.write("{}\t{}\n".format(et, count))
			fet.close()


			mList.append(thisM)
			mNames.append(et)
			
		#end if


	#end loop


	return mList, mNames

#end def ######## ######## ######## 
