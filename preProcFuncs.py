

import os.path
import sys




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
				keepGenes.append(lv[1])
			#end loop
		# Read in the edge types (# given in file)
		elif lv[0] == 'EDGE TYPES' :
			count = int(lv[2])
			for i in range(0, count) :
				line = f.readline()
				line = line.rstrip()
				lv = line.split('\t')
				keepEdges.append(lv[1])
			#end loop
		elif lv[0] == 'INDIRECT' :
			count = int(lv[2])
			for i in range(0, count) :
				line = f.readline()
				line = line.rstrip()
				lv = line.split('\t')
				keepEdges.append(lv[1])
			#end loop
		#end if
	#end loop

	return keepGenes, keepEdges, indirEdges

#end func ###################




