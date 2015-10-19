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


## initial length of buildAnArray
#addLen = 3000
#
## data type for edge array
#dt = np.dtype('a30')

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


# Initialize empty file
# This file stores what types were kept & how many edges
fn = path + oname + 'types.txt'
fet = open(fn, 'wb')
#fet.write("{}{}{}\n".format(et, delim, count))
fet.close()


####### ####### ####### ####### 
# DIRECT
# Create matrices from the gene-gene edges
for et in direct :
    if (et not in eTypes) :
        print "-- skipping: {}".format(et)
        continue
    #end if

    thisM = np.zeros([numG,numG])
    count = 0

#    for row in edgeArray :
#        if row[3] != et :
#            continue

    print "Creating matrix from edge type: {}".format(et)
    thisArray = edgeArray[edgeArray[:,3]==et]
    # increment entry at (i,j) = (gene0,gene1)
    for row in thisArray :
        thisM[geneDict[row[0]],geneDict[row[1]]] += 1
        thisM[geneDict[row[1]],geneDict[row[0]]] += 1
        count += 1
    #end loop

    # save to a file
    fn = path + oname + et
    print "    saving to {}".format(fn)
    np.save(fn, thisM)

#ERROR CHECK: save to a text file
#    fn = path + oname + et + '.txt'
#    print "    saving to {}".format(fn)
#    np.savetxt(fn, thisM, delimiter=delim)

    # This file stores what types were kept & how many edges
    fn = path + oname + 'types.txt'
    fet = open(fn, 'ab')
    fet.write("{}{}{}\n".format(et, delim, count))
    fet.close()
#end loop
del thisArray
del thisM




####### ####### ####### ####### 
# INDIRECT
# Create matrices from the node-gene edges
for et in indirect :
    if (et not in eTypes) :
        print "-- skipping: {}".format(et)
        continue
    #end if
    print "Building from: {}".format(et)

    # Create lists corresponding to the new edge types
    term_sm = list()
    term_md = list()
    term_lg = list()

    fn = path + ename + '.' + et + '.txt'
    fin = open(fn, 'rb')
    for line in fin :
        line = line.rstrip()
        lv = line.split(delim)

#        print lv[0], lv[1]

        # Determine the new edge type
        tSize = int(lv[1])
        if (tSize >= cutf[et][0]) & (tSize < cutf[et][1]) :
            term_sm.append(lv[0])
        elif  (tSize >= cutf[et][1]) & (tSize < cutf[et][2]) :
            term_md.append(lv[0])
        elif (tSize >= cutf[et][2]) & (tSize < cutf[et][3]) :
            term_lg.append(lv[0])
        #end if
    #end loop
    fin.close()


#    print term_sm

    # Create the first (small) matrix
    print "Creating matrix from edge type: {}".format(et+' < '+str(cutf[et][1]))
    thisM = np.zeros([numG,numG])
    count = 0
    for term in term_sm :
        # list of all edges with this term
        rowList = nodeDict[term]
        # Join all genes connected to term X
        for i in range(0, len(rowList)) :
            for j in range(0+1, len(rowList)) :
                # Find two genes joined by term
                # ASSUME: terms always in col 0, genes in col 1
                gA = edgeArray[i,1]
                gB = edgeArray[j,1]
                # Increment the entry(s) in the array (symmetric)
                thisM[geneDict[gA],geneDict[gB]] += 1
                thisM[geneDict[gB],geneDict[gA]] += 1
                count += 1
            #end loop
        #end loop
    #end loop
    if (count > 0) :
        # save to a file
        fn = path + oname + et + '_sm'
        print "    saving to {}".format(fn)
        np.save(fn, thisM)
        # This file stores what types were kept & how many edges
        fn = path + oname + 'types.txt'
        fet = open(fn, 'ab')
        fet.write("{}{}{}\n".format(et+'_sm', delim, count))
        fet.close()
    #end if

    # Create the second (medium) matrix
    print "Creating matrix from edge type: {}".format(et+' < '+str(cutf[et][2]))
    thisM = np.zeros([numG,numG])
    count = 0
    for term in term_md :
        # list of all edges with this term
        rowList = nodeDict[term]
        # Join all genes connected to term X
        for i in range(0, len(rowList)) :
            for j in range(0+1, len(rowList)) :
                # Find two genes joined by term
                # ASSUME: terms always in col 0, genes in col 1
                gA = edgeArray[i,1]
                gB = edgeArray[j,1]
                # Increment the entry(s) in the array (symmetric)
                thisM[geneDict[gA],geneDict[gB]] += 1
                thisM[geneDict[gB],geneDict[gA]] += 1
                count += 1
            #end loop
        #end loop
    #end loop
    if (count > 0) :
        # save to a file
        fn = path + oname + et + '_md'
        print "    saving to {}".format(fn)
        np.save(fn, thisM)
        # This file stores what types were kept & how many edges
        fn = path + oname + 'types.txt'
        fet = open(fn, 'ab')
        fet.write("{}{}{}\n".format(et+'_md', delim, count))
        fet.close()
    #end if

    # Create the second (medium) matrix
    print "Creating matrix from edge type: {}".format(et+' < '+str(cutf[et][3]))
    thisM = np.zeros([numG,numG])
    count = 0
    for term in term_lg :
        # list of all edges with this term
        rowList = nodeDict[term]
        # Join all genes connected to term X
        for i in range(0, len(rowList)) :
            for j in range(0+1, len(rowList)) :
                # Find two genes joined by term
                # ASSUME: terms always in col 0, genes in col 1
                gA = edgeArray[i,1]
                gB = edgeArray[j,1]
                # Increment the entry(s) in the array (symmetric)
                thisM[geneDict[gA],geneDict[gB]] += 1
                thisM[geneDict[gB],geneDict[gA]] += 1
                count += 1
            #end loop
        #end loop
    #end loop
    if (count > 0) :
        # save to a file
        fn = path + oname + et + '_lg'
        print "    saving to {}".format(fn)
        np.save(fn, thisM)
        # This file stores what types were kept & how many edges
        fn = path + oname + 'types.txt'
        fet = open(fn, 'ab')
        fet.write("{}{}{}\n".format(et+'_lg', delim, count))
        fet.close()
    #end if

#end loop


print "\nDone.\n"