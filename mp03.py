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



####### ####### ####### ####### 
# PARAMETERS


# INPUT ####
ename = 'hsa_dghmw_c'
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

print pTypes