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
# Output: <set_name>-<# of genes>.txt
#	sorted list of genes in the set
# ---------------------------------------------------------

import preProcFuncs as pp



####### ####### ####### ####### 
# PARAMETERS

# Name & location of file containing the samples
sFile = '_rawSets_enrichr.txt'
sDir = '../Dropbox/mp/samplesEnrichr'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION


# 1) Read the input file
if not sDir.endswith('/') :
	sDir = sDir + '/'
fname = sDir + sFile

# Store genes in a dict of lists
#	key: sampName; value: list(geneName)
sampDict = dict()
sampSet = set()

print("\nReading input file {}".format(sFile))
fin = open(fname, 'r')

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
fin.close()

# Get the unique list of sample names
sampNames = list(sampSet)
sampNames.sort()
sampSet.clear()
del sampSet



# 2) Save each set to a file
#		as a sorted list of gene names
print("Saving the sample files to {}".format(sDir))
for sn in sampNames :

	# Retrieve the list of genes in this sample
	geneList = sampDict[sn]
	geneList.sort()

	# Write the file
	fname = sDir + sn + '-{}.txt'.format(len(geneList))
	fout = open(fname, 'w')

	# each line contains a gene name, no empties at end
	lineCount = 0
	for gn in geneList :
		lineCount += 1
		if lineCount > 1 :
			fout.write('\n')
		fout.write(gn)
	#end loop
	fout.close()
#end loop
print("Wrote {} sample files.".format(len(sampNames)))



print("\nDone.\n")