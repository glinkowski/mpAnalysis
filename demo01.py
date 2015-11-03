# ----------------------------------------------------
# author: Greg Linkowski
# project: metapaths for KnowEnG
# 
# Metapths, demo part 01
#
# 
# ----------------------------------------------------
import numpy as np
import mpfuncs as mp



####### ####### ####### ####### 
# PARAMETERS

sname = ''
spath = '../samples/'

ename = 'toy_hsa_c'
epath = '../networks/'

mfile = ename + '.metapaths.txt'
gfile = ename + '.genes.txt'
rfile = ename + '.meta.gene_pairs.txt'
pfile = ename + '.meta.path_names.txt'
delim = '\t'

ofile = 'mp_demo01-'
opath = '../output/'

####### ####### ####### ####### 




####### ####### ####### ####### 
# BEGIN MAIN FUNCTION


# Collect the list of included genes
geneSet = set()
gf = open(path + gfile, 'rb')
for line in gf :
    lv = line.split('delim')
    geneSet.add(lv[0])
#end if
gf.close()


# Collect the list of included pairs
pairDict = dict()
row = 0
rf = open(path + rfile, 'rb')
for line in rf :
    lv = line.split('delim')
    pairDict[lv[0]] = row
    row += 1
#end if
#pairSet = set(pairDict.keys())
rf.close()


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
sGenes, sSize = mp.readWholeSample(path, sname, True, True)


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



# genes in the sample which aren't in the network
leftout = set()

# The main part of the show
sums = np.zeros([nCols])

numPairs = 0
for g1 in sGenes :
    if (g1 not in geneSet) :
        leftout.add(g1)
        continue
    #end if

    for g2 in sGenes :
        if (g2 not in geneSet) :
            leftout.add(g2)
            continue
        #end if

        if (g1 == g2) :
            continue
        #end if

        pair = g1 + '-' + g2
#        if (pair not in pairSet) :
#            continue
#        #end if

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