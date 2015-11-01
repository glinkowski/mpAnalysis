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

sname = ''

ename = 'toy_hsa_c'
path = '../networks/'


mfile = ename + '.metapath.txt'
gfile = ename + '.genes.txt'
pfile = ename + '.paths.txt'
delim = '\t'


indirect=(
    ['GO_term', 'motif_u5_gc', 'pfam_domain', 'KEGG'])

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


####### ####### ####### ####### 




####### ####### ####### ####### 
# BEGIN MAIN FUNCTION


## Collect the list of included genes
#geneSet = set()
#gf = open(path + gfile, 'rb')
#for line in gf :
#    lv = line.split('delim')
#    geneSet.add(lv[0])
##end if
#gf.close()


# Collect the list of included pairs
pairDict = dict()
row = 0
gf = open(path + gfile, 'rb')
for line in gf :
    lv = line.split('delim')
    pairDict[lv[0]] = row
    row += 1
#end if
pairSet = set(pairDict.keys())
gf.close()


# Collect the list of included paths
pTypes = list()
pf = open(path + pfile, 'rb')
for line in pf :
    lv = line.split('delim')
    pTypes.append(lv[0])
#end if
pf.close()


# Get the genes from the sample
print "Reading the sample {}".format(sname)
sGenes, sSize = mp.readWholeSample(path, sname)


# Load the metapaths matrix
print "Reading the list of metapths from {}".format(mfile)
nRows = len(pairDict)
nCols = sum( 1 for line in open(path + mfile, "rb") )
mpMatrix = np.zeros([nRows,nCols])

mf = open(path + mfile, 'rb')
col = 0
for line in mf :
    line = line.rstrip()
    lv = line.split('delim')

    row = 0
    for item in lv :
        mpMatrix[row, col] = item
        row += 1
    #end loop
    
    col += 1
#end loop




# The main part of the show
sums = np.zeros([nCols])

numPairs = 0
for g1 in sGenes :
#    if (g1 not in geneSet) :
#        continue
#    #end if

    for g2 in sGenes :
#        if (g2 not in geneSet) :
#            continue
#        #end if

        if (g1 == g2) :
            continue
        #end if

        pair = g1 + '-' + g2
        if (pair not in pairSet) :
            continue
        #end if

        numPairs += 1

        row = pairDict[pair]
        for col in range(0, nCols)
            sums[col] += mpMatrix[row, col]
        #end loop

        numPairs += 1
    #end loop
#end loop
avgs = sums / numPairs



# pull three random sets for comparison