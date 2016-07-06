# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# Pre-Processing of the enrichr sample files
#
# The enrichr samples come in a text file where all the
#	genes are jumbled. This script will extract each set
#	and output a file formatted like the MSigDB samples.
#
# Input: _rawSets_<source>.txt, tab-delimited
#	each line: set_name, gene_name, ?_value, source
# Output: <set_name>.txt
#	sorted list of genes in the set
# ---------------------------------------------------------

import preProcFuncs as pp
import time



####### ####### ####### ####### 
# PARAMETERS

# Name & location of file containing the samples
sFile = '_rawSets_enrichr.txt'
sDir = '../Dropbox/mp/samplesEnrichr'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()



# 1) Read the input file
if not sDir.endswith('/') :
	fname = sDir + '/' + sFile
else :
	fname = sDir + sFile
#end if

# store genes in a dict of lists
#	key: sampName; value: list(geneName)
sampDict = dict()
sampSet = set()

print("Reading input file {}".format(sFile))
fin = open(fname, 'r')
# #TODO: ?? write to gzip file (compress the file) -- maybe
# fzip = gzip.open(fname + '.gz', 'wb')

# Read each line in the file
lineCount = 0
for line in fin :
	lineCount += 1
	if len(line) < 1 :
		print("line {} is empty, exiting input file".format(lineCount))
	#end if

	# Extract the sample and gene names
	line = line.rstrip()
	lv = line.split('\t')
	sampName = lv[0]
	geneName = lv[1]

	# Save the gene into the appropriate list
	#	if list doesn't yet exist, create it
	if sampName in sampSet :
		sampDict[sampName].append(geneName)
	else :
		sampSet.add(sampName)
		sampDict[sampName] = [geneName]
#end loop

# Get the unique list of sample names
sampList = list(sampSet)
sampList.sort()
sampSet.clear()
del sampSet



