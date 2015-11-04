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

readMaxLines = 75
mpValDT = np.dtype(np.uint8)

numRandSets = 75
maxMPLen = 2

sname = 'POOLA_INVASIVE_BREAST_CANCER'
spath = '../samples/'

ofile = 'mp_demo01-g.txt'
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
geneDict = dict()
geneSet = set()
gf = open(epath + gfile, 'rb')
numG = 0
for line in gf :
    lv = line.split(delim)
    geneDict[lv[0]] = numG
    geneSet.add(lv[0])
    numG += 1
#end if
gf.close()


#- # Collect the list of included pairs
#- print "    reading the list of gene-pairs (column headers)"
#- pairDict = dict()
#- row = 0
#- rf = open(epath + rfile, 'rb')
#- for line in rf :
#-     lv = line.split(delim)
#-     pairDict[lv[0]] = row
#-     row += 1
#- #end if
#- #pairSet = set(pairDict.keys())
#- rf.close()


# Collect the list of included paths
print "    reading the list of calculated metapaths (row headers)"
pTypes = list()
pf = open(epath + pfile, 'rb')
for line in pf :
    lv = line.split(delim)
    pTypes.append(lv[0])
#end if
pf.close()


# Get the genes from the sample
print "Reading in the sample {}".format(sname)
sGenes, sSize = mp.readWholeSample(spath, sname, True, True)
sGenes.sort()


# Load the metapaths matrix
print "Reading the list of metapths from {}".format(mfile)
#nRows = len(pairDict)
nRows = numG * (numG - 1) / 2
#nCols = sum( 1 for line in open(path + mfile, "rb") )
nCols = len(pTypes)
#mpMatrix = np.zeros([nRows,nCols], dtype=np.dtype(np.uint16))
#mpMatrix = np.zeros([nRows,nCols], dtype='i16')
mpMatrix = np.zeros([nRows,nCols], dtype=mpValDT)

print "nRows = {}, nCols = {}".format(nRows, nCols)

mf = open(epath + mfile, 'rb')
col = 0
for line in mf :
    line = line.rstrip('\n')
    lv = line.split(delim)

#    print line
#    print lv[0], lv[1]
#    break

    row = 0
#    for item in lv :
    for i in range(0, len(lv)) :
#        mpMatrix[row, col] = int(item)
        mpMatrix[row, col] = int(float(lv[i]))
        row += 1
    #end loop
    
    col += 1

    print "    read line {}".format(col)

    if col == readMaxLines :
        nCols = col-1
        break

#    break
#end loop


#print mpMatrix[0,0], mpMatrix[0,5] #, mpMatrix[3,32]
#sys.exit()




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

        if ( len(leftout) == len(sGenes) ) :
            print "\nERROR: no sample genes in network!\n"
            sys.exit()
        #end if

#        if (g1 == g2) :
#            continue
#        #end if

#        pair = g1 + '-' + g2
#        if (pair not in pairSet) :
#            continue
#        #end if

        id1 = geneDict[g1]
        id2 = geneDict[g2]
        pairIndex = 0.5 * ( -pow(id1,2) + (id1*(2*numG + 1)) - (2*numG))
        pairIndex = pairIndex + (id2 - id1 - 1)


#        numPairs += 1

#        row = pairDict[pair]
        row = pairIndex
        for col in range(0, nCols) :
            sums[col] += mpMatrix[row, col]
        #end loop

        numPairs += 1
    #end loop
#end loop

#avgs = sums / numPairs
#avgs = np.divide(sums, numPairs)



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
    randList = list(randSet)
    del randSet
    randList.sort()


    for i in range(0, N) :
        g1 = randList[i]
        for j in range(i+1, N) :
            g2 = randList[j]

#            pair = g1 + '-' + g2
##TODO: this shouldn't be necessary ???
#            if (pair not in pairSet) :
#               continue
#            #end if

            id1 = geneDict[g1]
            id2 = geneDict[g2]
            pairIndex = 0.5 * ( -pow(id1,2) + (id1*(2*numG + 1)) - (2*numG))
            pairIndex = pairIndex + (id2 - id1 - 1)


            row = pairIndex
#            row = pairDict[pair]
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
#pairDict.clear()
del mpMatrix


# Write to the output file
of = open(opath + ofile, 'wb')

of.write("{}\n{}\n".format(sname, ename))
of.write("max MetaPath length: {}\n\n".format(maxMPLen))

of.write("  primary              {} random samples\n".format(numRandSets))
of.write("     sums        mean    st.d.     max       min     median\n")

for c in range(0, nCols) :
    rmean = np.mean(stats[c,:])
    rstd = np.std(stats[c,:])
    rmax = np.amax(stats[c,:])
    rmin = np.amin(stats[c,:])
    rmed = np.median(stats[c,:])

    of.write("{}\n".format(pTypes[c]))
    of.write("  {:>8.1f}   {:>8.1f} {:>7.1f} {:>9.1f} {:>9.1f} {:>9.1f}\n".format(
        sums[c], #avgs[c],
        rmean, rstd, rmax, rmin, rmed))
#end loop

# Make note of the left out genes
of.write("\nGenes in sample, but not in network:")
loList = list(leftout)
loList.sort()
for gene in loList :
    of.write("\n  {}".format(gene))
#end loop

of.close()