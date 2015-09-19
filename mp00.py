# ----------------------------------------------------
# author: Greg Linkowski
# project: metapaths for KnowEnG
# 
# Metapths, step 00
#
# Step through the edge file and fix known typos.
#   Save the corrections as a new file.
# Additionally, collect basic stats.
# ----------------------------------------------------




####### ####### ####### ####### 
# PARAMETERS

outname = 'all-v1'
path = '../../edgefiles'

infile = 'all_edges.txt'
outfile = outname + '.edge.txt'
delim = '\t'

####### ####### ####### ####### 




####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

typo1 = ['PPI_IntAct', 'PPI_intact']
correct1 = 'PPI_IntAct'

typo2 = ['STRING_cooccurrance', 'STRING_coocurrance']
correct2 = 'STRING_cooccurrence'

# satistics: Count # of each edge type
eTypes = dict()
# stats: count # of each node occurrence
nodes = dict()

ifile = open(path + infile, 'rb')
ofile = open(path + outfile, 'wb')

print "Checking {} for typos, saving to {}".format(ifile, ofile)

# Scan in each line, check for typos
#count = 0
for line in ifile :
#    count += 1

    # remove excess space & \n
    line = line.rstrip()
    lv = line.split(delim)

    # correct typo if it exists
    etype = lv[3]
    if etype in typo1 :
        etype = correct1
    elif etype in typo2 :
        etype = correct2
    #end if


    # stats: increment type count
    if etype in eTypes.keys() :
        eTypes[etype] += 1
    else :
        eTypes[etype] = 1
    #end if

    # stats: increment node count
    if lv[0] in nodes.keys() :
        nodes[lv[0]] += 1
    else :
        nodes[lv[0]] = 1
    #end if
    if lv[1] in nodes.keys() :
        nodes[lv[1]] += 1
    else :
        nodes[lv[1]] = 1
    #end if


    # write to output file
    ofile.write("{1}{0}{2}{0}{3}{0}{4}\n".format(delim,
        lv[0], lv[1], lv[2], etype))
#end loop
print "There were {} lines in the file.".format(count)

# close up shop
ifile.close()
ofile.close()




######## ######## ######## ########
# Calculate & output some basic graph stats
statfile = path + outname + '.stat.basic.txt'

print "Writing the stats file: {} ...".format(statfile)

sf = open(statfile, 'wb')


sf = open(outpath + outfile, "wb")
sf.write("---created by mp00.py : basic graph stats---\n\n")
sf.write("Basic edge & node statistics for the file:\n")
sf.write("\t{}\n\n".format(outfile))

# stats: edge counts
numEdges = sum(eTypes.values())
typeList = list(eTypes.keys())
typeList.sort()
numType = len(typeList)
sf.write("Total Edges: {:,}\n".format(numEdges))
for et in typeList :
    sf.write("  {:<20}{:>11,}{:>8.2f}%\n".format( et, etCount[et], 100 * etCount[et] / float(numEdges) ))
#end loop
sf.write("\n")

#stats: gene counts
gene_regex = ['ENSG', 'ENSM', 'FBgn', 'GB[45][0123456789]', 'WBGe', 'Y[ABCDEFGHIJKLMNOP][LR][0123456]', 'AT[12345CM]G']
gene_source = ['human', 'mouse', 'fly', 'bee', 'worm', 'yeast', 'grass']




# stats: node misc
nodeList = list(nodes.keys())
nodeList.sort()
sf.write("Total Nodes: {:,}\n".format( len(NodeSet) ))
sf.write("Node Degree distributions by type:\n")
sf.write("    -  Node Type  -      -# Nodes-     -Mean-  \n")
#of.write("         -  Node Type  -          -# Nodes-     -Mean- -StdDev-  -Min-  -Max- \n")
# first the genes
#of.write(" {:<11} {:>8,} {:9.2f} {:8.2f} {:>5} {:>7}\n".format("genes", len(stats), np.mean(stats), np.std(stats), np.min(stats), np.max(stats) ))


#numGenes = len(np.unique(Edges[:,1]))
#of.write(" {:<20} {:>11,} {:>11.1f} \n".format("genes", numGenes, numEdges / float(numGenes) ))
# then the rest
for et in typeList :
    numType = len(np.unique(Edges[Edges[:,3] == et][:,0]))
    sf.write(" {:<20} {:>11,} {:>11.1f} \n".format( et, numType, numEdges / float(numType)))
#end loop
sf.write("\n")

# stats: go terms
if "GO_term" in typeList :
    # Extract the go terms
    goList = nodeList[nodeList[0:3] == "GO:"]
    numgo = len(goList)

    # Get the size of each go term
    goLen = dict()
    goSet = set()
    for gt in goList :
        if gt in goSet :
            goLen[gt] += 1
        else :
            goSet.add(gt)
            goLen[gt] = 1
        #end if
    #end loop

    # Save these to a file
    gofile = path + outname + ".goList.txt"
    print "Saving GO terms to {}".format(gofile)
    gf = open(gofile, 'wb')
    goList = list(goSet)
    goList.sort()
    for gt in goList :
        gf.write("{}{}{}\n".format(gt, delim, goLen[gt]))
    gf.close()
    del goSet

    sf.write("Number of GO terms of length ...")
    # Count how many terms exist of each length
    lenList = list(goLen.values())
    for l in range(0,11) :
        a = l * 10
        b = a + 9
        gcount = ((lenList >= a) & (lenList >= b)).sum()
        sf.write("  {:>3} -{:>3} :  {:>4}\n".format(a, b, gcount))
    #end loop
    for l in range(2,7) :
        a = l * 50
        b = a + 49
        gcount = ((lenList >= a) & (lenList >= b)).sum()
        sf.write("  {:>3} -{:>3} :  {:>4}\n".format(a, b, gcount))
    #end loop
    sf.write("  300+     :  {:>4}\n".format(a, b, gcount))
    sf.write('\n')

    del goList
    goLen.clear()
#end if

sf.close()




print "\nDone.\n"
####### ####### ####### ####### 