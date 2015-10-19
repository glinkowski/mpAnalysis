# ----------------------------------------------------
# author: Greg Linkowski
# project: metapaths for KnowEnG
# 
# Metapths, step 02
#
# Create a network consisting only of gene-gene edges.
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
outfile = ename + '.edge_gonly.txt'
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
cutf['GO_term'] = [0, 50, 100, top]
cutf['motif_u5_gc'] = [0, 3000, 8000, top]
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

# Create a new array for Indirect edge types

fout = open(path+outfile, 'wb')
#buildAnArray = np.empty((addLen,4), dtype=dt)
index = 0
#delLen = 0
for etx in range(0,4) :
    cType = indirect[0][etx]
#for et in indirect[1,:] :
    # only altering the Indirect edges
    if cType not in eTypes :
#    if et not in eTypes :
        continue
    #end if

    # Load the term dictionary from the file
    termDict = dict()
    fn = ename + '.' + indirect[1][etx] + '.txt'
    print "Reading term file {}".format(fn)
    fdic = open(path+fn, 'rb')
    for line in fdic :
        line = line.rstrip()
        lv = line.split(delim)
        termDict[lv[0]] = lv[1]
    #end if
    fdic.close()


#    # Extract the edges of type X
#    tempArray = edgeArray[edgeArray[:,3] == indirect[0,etx]]

    print "  defining new edges for {}".format(cType)
    # Get a list of all terms of type X
    termList = list(termDict.keys())

    for term in termList :
        # get the indices where term X appears
        rowList = nodeDict[term]
#        delLen += len(rowList)

        # Determine the new edge type
        tSize = termDict[term]
        if (tSize >= cutf[cType][0]) & (tSize < cutf[cType][1]) :
            newEdge = cType + '_sm'
        elif  (tSize >= cutf[cType][1]) & (tSize < cutf[cType][2]) :
            newEdge = cType + '_md'
        elif (tSize >= cutf[cType][2]) & (tSize < cutf[cType][3]) :
            newEdge = cType + '_lg'
        #end if

        # Join all genes connected to term X
        for i in range(0, len(rowList)) :
            for j in range(0+1, len(rowList)) :

#                # Adjust array size if needed
#                if (index >= addLen) :
##                if (index >= size(buildAnArray)[0]) :
#                    newArray = np.empty((addLen, 4), dtype=dt)
#                    buildAnArray = np.vstack( (buildAnArray, newArray) )
#                    del newArray
##                    arrLen = 2 * arrLen + 1
#                #end if

                # Add new gene-gene edge to network file
                if index == 0 :
                    fout.write("{1}{0}{2}{0}{3}{0}{4}".format(delim, 
                        edgeArray[rowList[i], 1], 
                        edgeArray[rowList[j], 1], 
                        min(edgeArray[rowList[i]]), newEdge))
                else :
                    fout.write("\n{1}{0}{2}{0}{3}{0}{4}".format(delim, 
                        edgeArray[rowList[i], 1], 
                        edgeArray[rowList[j], 1], 
                        min(edgeArray[rowList[i]]), newEdge))
                #end if


 #               # Add new gene-gene edge to array
 #               # Assumption: genes will always be in col 2
##                print edgeArray[rowList[i],0]
##                buildAnArray[index][0] = 'oats'
 #               buildAnArray[index,0] = edgeArray[rowList[i], 1] 
 #               buildAnArray[index,1] = edgeArray[rowList[j], 1]
 #               buildAnArray[index,2] = str(min(edgeArray[rowList[i]][2]))
#                buildAnArray[index,3] = newEdge
##                buildAnArray[index] = ([
##                    edgeArray[rowList[i],1], 
##                    edgeArray[rowList[j],1], 
##                    str(min(edgeArray[rowList[i],2])), 
##                    newEdge])
                index += 1
            #end loop 'j'
        #end loop 'i'
    #end loop 'term'

    termDict.clear()
#end loop 'etx'

dirSet = set(direct)
for row in edgeArray :
    fout.write("\n{1}{0}{2}{0}{3}{0}{4}".format(delim, 
    row[0], row[1], row[2], row[3]))
#end if

fout.close()

#print "Finished creating new edges ..."
## Delete the unused rows from buildAnArray
#arrLen = np.size(buildAnArray)[0]
#buildAnArray = np.delete(buildAnArray, np.s_[index:arrLen],0)
#
#print "  combining into final network"
## Efficient way to do this?
## Build final array from desired rows, and buildAnArray
#keepRows = np.size(edgeArray)[0] - delLen
#finArray = np.empty(( (index + keepRows), 4))
#dirSet = set(direct)
#r = 0
##for r in range(0, keepRows) :
#for row in edgeArray :
#    if row[3] in dirSet :
#        finArray[r] = row
#        r += 1
#    #end if
##end loop
#del edgeArray
#for row in buildAnArray :
#    finArray[r] = row
#    r += 1
##end loop
#del buildAnArray
#
#
## TODO:
## Build new node dictionary
#
#
#
## Output to file
#print "Writing new edge file to {}".format(outfile)
#mp.writeEdgeFile(finArray, path+outfile, delim)


print "\nDone.\n"