# ----------------------------------------------------
# author: Greg Linkowski
# project: metapaths for KnowEnG
# 
# pre-Metapths, create smaller toy network
#
# Goal: Create a smaller output file for demos
# 	Remove XX% of the nodes in the original network
# ----------------------------------------------------
import sys
#import re
import mpfuncs as mp
import random



####### ####### ####### ####### 
# PARAMETERS

ename = 'hsa_dghmw_c'
#ename = 'all-v1'
oname = 'toy2_hsa'
path = '../networks/'

infile = ename + '.edge.txt'
outfile = oname + '.edge.txt'
outdict = oname + '.dict.txt'
delim = '\t'

# percent of nodes to remove
prune = 50

## species               gene designation
#include_human = True   # ENSG
#include_mouse = True   # ENSM
#include_fly = True     # FBgn
#include_bee = True     # GB[45][0123456789]
#include_worm = True    # WBGe
#include_yeast = True   # Y[ABCDEFGHIJKLMNOP][LR][0123456]
#include_grass = True   # AT[12345CM]G

######## ######## ######## ########



######### ######## ######## ########
## IMPLEMENT PARAMETERS
#
#gene_regex = ['ENSG', 'ENSM', 'FBgn', 'GB[45][0123456789]', 'WBGe', 'Y[ABCDEFGHIJKLMNOP][LR][0123456]', 'AT[12345CM]G']
#gene_select = [include_human, include_mouse, include_fly, include_bee, include_worm, include_yeast, include_grass]
#
#keepGenes = list()
#for i in range(0, len(gene_regex)) :
#    if gene_select[i] :
#        keepGenes.append(gene_regex[i])
##end loop
#
##print keepGenes
#
######### ######## ######## ########



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

if (path+outfile) == (path+infile) :
    print "ERROR: Input and output are same file: {}".format(path+outfile)
    sys.exit()
#end if


print "Reading in the original edge file..."
edgeArray, nodeDict = mp.readEdgeFile(path+infile, delim)


# Build the list of nodes to exclude
remSet = set()
# normalize prune value to [0,1]
check = prune / 100.0

nodeList = list(nodeDict.keys())
for node in nodeList :
	# randomly add nodes to remove list
	test = random.random()
	if test <= check :
		remSet.add(node)
	#end if
#end loop

print "Saving the reduced edge file as {}".format(outfile)
f = open(path + outfile, 'wb')
firstLine = True
for edge in edgeArray :

	# Skip nodes in the remSet
	if edge[0] in remSet :
		continue
	elif edge[1] in remSet :
		continue
	#end if

	# Save the rest of the edges to the new file
	if firstLine :
		f.write("{1}{0}{2}{0}{3}{0}{4}".format(delim, edge[0], edge[1], edge[2], edge[3]))
		firstLine = False
	else :
		f.write("\n{1}{0}{2}{0}{3}{0}{4}".format(delim, edge[0], edge[1], edge[2], edge[3]))
	#end if
#end loop
f.close()


print "\nDone.\n"