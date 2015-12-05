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

import sys


####### ####### ####### ####### 
# PARAMETERS

numRandSets = 3
maxMPLen = 2

sname = 'LEE_LIVER_CANCER_E2F1'
spath = '../samples/'

ofile = 'mp_demo01-a.txt'
opath = '../output/'

ename = 'toy_hsa_c'
epath = '../networks/'

mfile = ename + '.metapaths.txt'
gfile = ename + '.genes.txt'
rfile = ename + '.meta.gene_pairs.txt'
pfile = ename + '.meta.path_names.txt'
delim = '\t'

####### ####### ####### ####### 




####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

print "Loading files ..."

# Collect the list of included genes
print "    loading genes in network"
geneSet = set()
gf = open(epath + gfile, 'rb')
for line in gf :
    lv = line.split('delim')
    geneSet.add(lv[0])
#end if
gf.close()


# Collect the list of included pairs
print "    reading the list of gene-pairs (column headers)"
pairDict = dict()
row = 0
rf = open(epath + rfile, 'rb')
for line in rf :
    lv = line.split('delim')
    pairDict[lv[0]] = row
    row += 1
#end if
#pairSet = set(pairDict.keys())
rf.close()


# Collect the list of included paths
print "    reading the list of calculated metapaths (row headers)"
pTypes = list()
pf = open(epath + pfile, 'rb')
for line in pf :
    lv = line.split('delim')
    pTypes.append(lv[0])
#end if
pf.close()


# Get the genes from the sample
print "Reading in the sample {}".format(sname)
sGenes, sSize = mp.readWholeSample(spath, sname, True, True)
sGenes.sort()


# Load the metapaths matrix
print "Reading the list of metapths from {}".format(mfile)
nRows = len(pairDict)
#nCols = sum( 1 for line in open(path + mfile, "rb") )
nCols = len(pTypes)
mpMatrix = np.zeros([nRows,nCols], dtype=np.uint16)

print "nRows = {}, nCols = {}".format(nRows, nCols)

mf = open(epath + mfile, 'rb')
col = 0
for line in mf :
    line = line.rstrip('\n')
    lv = line.split('delim')

#    print line
#    print lv[0], lv[1]
#    break

    row = 0
#    for item in lv :
    for i in range(0, len(lv)) :
#        mpMatrix[row, col] = int(item)
        mpMatrix[row, col] = int(lv[i])
        row += 1
    #end loop
    
    col += 1
#end loop

print mpMatrix[0,0], mpMatrix[0,5], mpMatrix[3,32]



sys.exit()




####### ####### ####### ####### 
# PRIMARY SAMPLE

# genes in the sample which aren't in the network
leftout = set()

# The main part of the show
sums = np.zeros([nCols])

print "Searching for metapaths in the sample"
numPairs = 0
#for g1 in sGenes :
for i in range(0, len(sGenes)) :
    g1 = sGenes[i]
    if (g1 not in geneSet) :
        leftout.add(g1)
        continue
    #end if

#    for g2 in sGenes :
    for j in range(i+1, len(sGenes)) :
        g2 = sGenes[j]
        if (g2 not in geneSet) :
            leftout.add(g2)
            continue
        #end if

#        if (g1 == g2) :
#            continue
#        #end if

        pair = g1 + '-' + g2
#        if (pair not in pairSet) :
#            continue
#        #end if

        numPairs += 1

        row = pairDict[pair]
        for col in range(0, nCols) :
            sums[col] += mpMatrix[row, col]
        #end loop

        numPairs += 1
    #end loop
#end loop
avgs = sums / numPairs



####### ####### ####### ####### 
# RANDOM SAMPLES

# pull XX random sets for comparison

# Create array to hold stats
stats = np.zeros([nCols, numRandSets])

print "Searching for metapaths in {} random samples".format(numRandSets)
for r in range(0, numRandSets) :

    # Get a random set of genes
    N = len(sGenes) - len(leftout)
    randSet =  mp.selectRandomNodes(N, geneSet)


    for i in range(0, len(sGenes)) :
        g1 = sGenes[i]
        for j in range(i+1, len(sGenes)) :
            g2 = sGenes[j]

            pair = g1 + '-' + g2
#TODO: this shouldn't be necessary ???
            if (pair not in pairSet) :
               continue
            #end if

            row = pairDict[pair]
            for col in range(0, nCols) :
                stats[col, r] += mpMatrix[row, col]
            #end loop
        #end loop
    #end loop
#end loop



####### ####### ####### ####### 
# OUTPUT

print "Done searching. Writing output file."

# Free some memory
pairDict.clear()
del mpMatrix


# Write to the output file
of = open(opath + ofile, 'wb')

of.write("{}\n{}\n".format(sname, ename))
of.write("max MetaPath lenght: {}}\n\n".format(maxMPLen))

of.write("  primary      {} random samples\n")
of.write("    mean       mean    st.d.     max       min     median\n")

for c in range(0, nCols) :
    rmean = np.mean(stats[c,:])
    rstd = np.std(stats[c,:])
    rmax = np.amax(stats[c,:])
    rmin = np.amin(stats[c,:])
    rmed = np.median(stats[c,:])

    of.write("{}\n".format(pTypes[c]))
    of.write("  {:>8}   {:>8} {:>7} {:>9} {:>9} {:>9}\n".format(
        avgs[c], rmean, rstd, rmax, rmin, rmed))
#end loop
of.close()