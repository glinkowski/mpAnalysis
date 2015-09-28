# ----------------------------------------------------
# author: Greg Linkowski
# project: metapaths for KnowEnG
# 
# Metapths, step 01
#
# Normalize weights per edge type.
# ----------------------------------------------------

import mpfuncs.py as mp
import numpy as np



####### ####### ####### ####### 
# PARAMETERS

ename = 'all-v1'
path = '../networks/'

infile = ename + '.edge.txt'
outfile = ename + '.edge_norm.txt'
#outfile = ename + '.gene-only.txt'
delim = '\t'


#indirect=['GO_term', 'motif_u5_gc', 'pfam_domain'] 
#direct = ['prot_homol']

####### ####### ####### ####### 




####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

print "Reading in the non-normalized edge file..."
edgeArray, nodeDict =  mp.readEdgeFile(path+infile, delim)

# get the unique edge types in the network
eTypes = list( np.unique(edgeArray[:,3]) )
eTypes.sort()

print "{} contains the following edges:".format(infile)
print eTypes

for et in eTypes :

    # not sure this line works
    # may need to make temp array
    values = edgeArray[edgeArray[:,3]==et, 3]

    vmin = np.amin(values)
    vmax = np.amin(values)
    vdif = vmax - vmin

    for i in len(edgeArray)[0] :
        if edgeArray[i,3] == et :

            if (vdif == 0) :
                # if vmin=vmax, then there's only one value
                temp = 1
            else :
                temp = float((edgeArray[i,2] - vmin) / vdif)
            #end if

            edgeArray[i,2] = str(temp)
        #end if
    #end loop

#end loop

print "Writing the normalized edge file to {}".format(outfile)
mp.writeEdgeFile(edgeArray, path+outfile, delim)



print "\nDone.\n"