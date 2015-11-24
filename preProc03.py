# ----------------------------------------------------
# author: Greg Linkowski
# project: metapaths for KnowEnG
# 
# Metapths, step 03
#
# From the matricies, caluculate metapaths.
#	The list of paths for this network is stored in
#	file: ...matrix.types.txt
#	Step through the combinations in this file to
#	calculate & store the output.
# ----------------------------------------------------
import numpy as np
import mpfuncs as mp

import sys


####### ####### ####### ####### 
# PARAMETERS


# INPUT ####
ename = 'toy_hsa_c'
#ename = 'hsa_dghmw_c'
#ename = 'all-v1'
path = '../networks/'
delim = '\t'

# The matrices for each edge type
fmatrix = ename + '.matrix.'
# The row/column names -- list of genes
fgenes = ename + '.genes.txt'
# The types of matrices made from this network
ftypes = ename + '.matrix.types.txt'


# OUTPUT ####
# The output matrix -- metapath counts
omatrix = ename + '.metapaths.'
# output column names -- the gene pairs
ogpairs = ename + '.meta.gene_pairs.txt'
# output row names -- the meta paths
opnames = ename + '.meta.path_names.txt'



#indirect=(
#    ['GO_term', 'motif_u5_gc', 'pfam_domain', 'KEGG'])
#
#direct = (['prot_homol', 'STRING_experimental',
#    'STRING_coexpression',  'STRING_textmining', 
#    'STRING_neighborhood', 'PPI_IntAct', 
#    'STRING_database', 'STRING_cooccurrence', 
#    'STRING_fusion', 'PPI_MINT', 'PPI_DIP', 'PPI_BioGRID'])

#
#top = 65534
#cutf = dict()
#cutf['GO_term'] = [0, 40, 100, 450]
#cutf['motif_u5_gc'] = [0, 3000, 8000, 10000]
#cutf['pfam_domain'] = [0, 40, 300, top]
#cutf['KEGG'] = [0, 20, 100, top]


####### ####### ####### ####### 




####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

# Read in the paths collected from network
pTypes = list()
ftemp = open(path + ftypes, 'rb')
for line in ftemp :
	line = line.rstrip()
	lv = line.split(delim)
	pTypes.append(lv[0])
#end loop


# Read in the list of genes in this network
#   dictionary will indicate what index to use in matrix
print "Reading in the list of genes..."
geneDict = dict()
gf = open(path + fgenes, 'rb')
numG = 0
for line in gf :
    lv = line.split(delim)
    geneDict[lv[0]] = numG
    numG += 1
#end loop
gf.close()


print "Types included in this network file: "
print " ", pTypes






# OUTPUT: list of gene pairs in order --> column names
geneDict.clear() # no longer needed
# Get the ordered list of genes in these matrices
gList = list()
gf = open(path + fgenes, 'rb')
for line in gf :
    lv = line.split(delim)
    gList.append(lv[0])
#end loop
gf.close()

print gList[0], gList[30], gList[8]

# Build the ordered list of gene-pairs
gpList = list()
for i in range(0, numG) :
	for j in range(i+1, numG) :
		gpList.append(gList[i] + '-' + gList[j])
	#end loop
#end loop

print gpList[1], gpList[30]

# Save the list of gene-pairs to a file
rfile = open(path + ogpairs, 'wb')
firstLine = True
for item in gpList :
	if (not firstLine) :
		rfile.write("\n")
	#end if
	firstLine = False
	rfile.write("{}".format(item))
#end loop
rfile.close()



#print "\nExiting early\n"
#sys.exit()





# Preparation for saving metapaths
print "Preparing to calculate metapaths ..."
mPaths = list()
mpfn = path + omatrix + 'txt'
mpfile = open(mpfn, 'wb')
mpfile.close()

# Create the length-1 metapath list
firstLine = True
for mp in pTypes :

	print "    finding {}".format(mp)
	# Load the gene-gene matrix from file
	A = np.load(path + fmatrix + mp + '.npy')
	L = A.shape[0]

	# the row of values to write to master mp file
	row = list()

#	count = 0
	for i in range(0, L) :
		for j in range(i+1, L) :

#			count += A[i,j]
			row.append(A[i,j])
		#end loop
	#end loop


	print "        saving..."
	# OUTPUT: write the row to the file
	mpfile = open(mpfn, 'ab')
	if (not firstLine) :
		mpfile.write("\n")
	#end if
	firstLine = False

	firstItem = True
	for r in row :
		if (not firstItem) :
			mpfile.write("{}".format(delim))
		#end if
		mpfile.write("{}".format(r))
		firstItem = False
	#end loop
	mpfile.close()


	# track: the name of this metapath & count (in order)
	mPaths.append([mp, sum(row)])
#end loop


# Create the length-2 metapath list
for mp1 in pTypes :
	for mp2 in pTypes :
		mp = mp1 + '-' + mp2

		print "    finding {}".format(mp)
		# Load the gene-gene matrices from file
		A = np.load(path + fmatrix + mp1 + '.npy')
		L = A.shape[0]
		B = np.load(path + fmatrix + mp2 + '.npy')

		# Calculte the new matrix
		C = np.multiply(A, B)

		# the row of values to write to master mp file
		row = list()
		for i in range(0, L) :
			for j in range(i+1, L) :
				row.append(C[i,j])
			#end loop
		#end loop

		print "        saving..."
		# OUTPUT: write the row to the file
		mpfile = open(mpfn, 'ab')
		mpfile.write("\n")

		firstItem = True
		for r in row :
			if (not firstItem) :
				mpfile.write("{}".format(delim))
			#end if
			mpfile.write("{}".format(int(r)))
			firstItem = False
		#end loop
		mpfile.close()

		# track: the name of this metapath & count (in order)
		mPaths.append([mp, sum(row)])
	#end loop
#end loop





## OUTPUT: list of gene pairs in order --> column names
#geneDict.clear() # no longer needed
## Get the ordered list of genes in these matrices
#gList = list()
#gf = open(path + fgenes, 'rb')
#for line in gf :
#    lv = line.split(delim)
#    gList.append(lv[0])
##end loop
#gf.close()
#
## Build the ordered list of gene-pairs
#gpList = list()
#for i in range(0, numG) :
#	for j in range(i+1, numG) :
#		gpList.append(gList[i] + '-' + gList[j])
#	#end loop
##end loop
#
## Save the list of gene-pairs to a file
#rfile = open(path + ogpairs, 'wb')
#firstLine = True
#for item in gpList :
#	if (not firstLine) :
#		pfile.write("\n")
#	#end if
#	firstLine = False
#	rfile.write("{}{}{}".format(item[0]))
##end loop
#rfile.close()



# OUTPUT: list of metapaths in order --> row names
pfile = open(path + opnames, 'wb')
firstLine = True
for item in mPaths :
	if (not firstLine) :
		pfile.write("\n")
	#end if
	firstLine = False

	pfile.write("{}{}{}".format(item[0],delim,item[1]))
#end loop
pfile.close()


print "\nDone. ... for now\n"