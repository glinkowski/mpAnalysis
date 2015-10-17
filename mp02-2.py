# ----------------------------------------------------
# author: Greg Linkowski
# project: metapaths for KnowEnG
# 
# Metapths, step 02
#
# Create matrices consisting only of gene-gene edges.
#   Replace indirect edges with direct edges. Along
#   the way, create new edges reflective of the
#   strength of membership for that term.
# For example, 3-4 new GO_term edge types to represent
#   higher import of sharing a smaller term.
# ----------------------------------------------------
import numpy as np
import mpfuncs as mp



####### ####### ####### ####### 
# PARAMETERS

ename = 'hsa_dghmw_c'
#ename = 'all-v1'
path = '../networks/'
infile = ename + '.edge_norm.txt'
#outfile = ename + '.edge_gonly.txt'
oname = ename + '.matrix.'
delim = '\t'


indirect=(
    ['GO_term', 'motif_u5_gc', 'pfam_domain', 'KEGG'], 
    ['goList', 'mtList', 'pfList', 'kgList']
    )

direct = (['prot_homol', 'STRING_experimental',
    'STRING_coexpression',  'STRING_textmining', 
    'STRING_neighborhood', 'PPI_IntAct', 
    'STRING_database', 'STRING_cooccurrence', 
    'STRING_fusion', 'PPI_MINT', 'PPI_DIP', 'PPI_BioGRID'])


top = 65534
cutf = dict()
cutf['GO_term'] = [0, 40, 100, 450]
cutf['motif_u5_gc'] = [0, 3000, 8000, 10000]
cutf['pfam_domain'] = [0, 40, 300, top]
cutf['KEGG'] = [0, 20, 100, top]


# initial length of buildAnArray
addLen = 3000

# data type for edge array
dt = np.dtype('a30')

####### ####### ####### ####### 




####### ####### ####### ####### 
# BEGIN MAIN FUNCTION


# Load the network from file
print "Reading in the edge file..."
edgeArray, nodeDict =  mp.readEdgeFile(path+infile, delim)

# Get the unique edge types in the network
eTypes = list( np.unique(edgeArray[:,3]) )
eTypes.sort()
print "{} contains the following edges:".format(infile)
print eTypes

# Get the list of genes in this network
#   dictionary will indicate what index to use in matrix
#ASSUMPTION: file is sorted
print "Reading in the list of genes..."
geneDict = dict()
gf = open(path + ename + '.genes.txt', 'rb')
numG = 0
for line in gf :
    lv = line.split(delim)
    geneDict[lv[0]] = numG
    numG += 1
#end loop
gf.close()


# Create matrices from the gene-gene edges
for et in direct :
    if (et not in eTypes) :
        continue
    #end if

    thisM = np.zeros([numG,numG])

#    for row in edgeArray :
#        if row[3] != et :
#            continue


    thisArray = edgeArray[edgeArray[:,3]==et]
    # increment entry at (i,j) = (gene0,gene1)
    for row in thisArray :
        thisM[geneDict[row[0]],geneDict[row[1]]] += 1
        thisM[geneDict[row[1]],geneDict[row[0]]] += 1
    #end loop

    # save to a file
    fn = path + oname + et + '.txt'
    np.save(fn, thisM)

#ERROR CHECK: save to a text file
    np.savetxt(fn, thisM, delimiter=delim)

#end loop
del thisM


# Create matrices from the gene-gene edges
for et in indirect :
    if (et not in eTypes) :
        continue
    #end if

#TODO: create these matrices


#end loop