# ----------------------------------------------------
# author: Greg Linkowski
# project: metapaths for KnowEnG
# 
# Metapths, step 01
#
# Normalize weights per edge type.
# ----------------------------------------------------

import mpfuncs as mp
import numpy as np



####### ####### ####### ####### 
# PARAMETERS

ename = 'hsa_dghmw_c'
#ename = 'all-v1'
path = '../networks/'

infile = ename + '.edge.txt'
outfile = ename + '.edge_norm.txt'
outdict = ename + '.dict.txt'
delim = '\t'

####### ####### ####### ####### 




####### ####### ####### ####### 
# BEGIN MAIN FUNCTION


print "Reading in the non-normalized edge file..."
edgeArray, nodeDict = mp.readEdgeFile(path+infile, delim)


# get the unique edge types in the network
eTypes = list( np.unique(edgeArray[:,3]) )
eTypes.sort()
print "{} contains the following edges:".format(infile)
print eTypes


print "Normalizing edge type ..."
for et in eTypes :
    print "    {}".format(et)

    # extract the edge weights for only this edge type
    values = edgeArray[edgeArray[:,3]==et, 2]
    # convert entries from str to float
    v2 = np.zeros((len(values),1))
    for i in range(0,len(values)) :
        v2[i] = float(values[i])
    #end loop

    vmin = float(np.amin(v2))
    vmax = float(np.amax(v2))
    vdif = vmax - vmin
#    print vmin, vmax, vdif

    for i in range(0, edgeArray.shape[0]) :
        # calculate normalized value & update array
        if edgeArray[i,3] == et :
            if (vdif == 0) :
                # if vmin=vmax, then there's only one value
                temp = 1
            else :
                # else change value to range [0,1]
                temp = (float(edgeArray[i,2]) - vmin) / vdif
            #end if
            edgeArray[i,2] = str(temp)
        #end if
    #end loop

#end loop


print "Writing the normalized edge file to {}".format(outfile)
mp.writeEdgeFilePlus(edgeArray, nodeDict, path+outfile, path+outdict, delim)



print "\nDone.\n"