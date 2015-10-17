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
import sys
import re



####### ####### ####### ####### 
# PARAMETERS

#outname = 'hsa_dghmw_c'
outname = 'all-v1'
path = '../networks/'

#infile = 'hsa_dghmw.edge.txt'
infile = 'all_edges.txt'
outfile = outname + '.edge.txt'
delim = '\t'

# species               gene designation
include_human = True   # ENSG
include_mouse = True   # ENSM
include_fly = True     # FBgn
include_bee = True     # GB[45][0123456789]
include_worm = True    # WBGe
include_yeast = True   # Y[ABCDEFGHIJKLMNOP][LR][0123456]
include_grass = True   # AT[12345CM]G

######## ######## ######## ########



######## ######## ######## ########
# IMPLEMENT PARAMETERS

gene_regex = ['ENSG', 'ENSM', 'FBgn', 'GB[45][0123456789]', 'WBGe', 'Y[ABCDEFGHIJKLMNOP][LR][0123456]', 'AT[12345CM]G']
gene_select = [include_human, include_mouse, include_fly, include_bee, include_worm, include_yeast, include_grass]

keepGenes = list()
for i in range(0, len(gene_regex)) :
    if gene_select[i] :
        keepGenes.append(gene_regex[i])
#end loop

#print keepGenes

######## ######## ######## ########



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

if (path+outfile) == (path+infile) :
    print "ERROR: Input and output are same file: {}".format(path+outfile)
    sys.exit()
#end if

typo1 = ['PPI_IntAct', 'PPI_intact']
correct1 = 'PPI_IntAct'

typo2 = ['STRING_cooccurrance', 'STRING_coocurrance']
correct2 = 'STRING_cooccurrence'

# satistics: Count # of each edge type
eTypes = dict()
etSet = set()
# stats: count # of each node occurrence
nodes = dict()
nodeSet = set()
# stats: get size of each GO term
goTerms = dict()
goSet = set()
# stats: get size of each MT term
mtTerms = dict()
mtSet = set()
# stats: get size of each Protein Family
pfTerms = dict()
pfSet = set()
# stats: get size of each KEGG pathway
kgTerms = dict()
kgSet = set()
# stats: get degree for each gene
geneDict = dict()
geneSet = set()


ifile = open(path + infile, 'rb')
ofile = open(path + outfile, 'wb')

print "Checking {} for typos, saving to {}".format(infile, outfile)

# Scan in each line, check for typos
count = 0
for line in ifile :
    if count > 0 :
        ofile.write('\n')
    #end if

    count += 1
    if (count % 10000000) == 0 :
        print "    on line {:,}".format(count)


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
    if etype in etSet :
        eTypes[etype] += 1
    else :
        etSet.add(etype)
        eTypes[etype] = 1
    #end if

    # stats: increment node count
    if lv[0] in nodeSet :
        nodes[lv[0]] += 1
    else :
        nodeSet.add(lv[0])
        nodes[lv[0]] = 1
    #end if
    if lv[1] in nodeSet :
        nodes[lv[1]] += 1
    else :
        nodeSet.add(lv[1])
        nodes[lv[1]] = 1
    #end if

    # stats: increment GO Term count
    if etype == 'GO_term' :
        # which node is the go term?
        if lv[0][0:3] == 'GO:' :
            term = lv[0]
        else :
            term = lv[1]
        #end if
        # increment that term's occurrence
        if term in goSet :
            goTerms[term] += 1
        else :
            goTerms[term] = 1
            goSet.add(term)
        #end if
    #end if

    # stats: increment MT term count
    if etype == 'motif_u5_gc' :
        # which node is the go term?
        if lv[0][0:3] == 'MT_' :
            term = lv[0]
        else :
            term = lv[1]
        #end if
        # increment that term's occurrence
        if term in mtSet :
            mtTerms[term] += 1
        else :
            mtTerms[term] = 1
            mtSet.add(term)
        #end if
    #end if

    # stats: increment P Fam count
    if etype == 'pfam_domain' :
        # which node is the go term?
        if lv[0][0:5] == 'Pdom_' :
            term = lv[0]
        else :
            term = lv[1]
        #end if
        # increment that term's occurrence
        if term in pfSet :
            pfTerms[term] += 1
        else :
            pfTerms[term] = 1
            pfSet.add(term)
        #end if
    #end if

    # stats: increment KEGG pathway count
    if etype == 'KEGG' :
#        # which node is the KEGG term?
#        if lv[0][0:3] == 'kg:' :
#            term = lv[0]
#        else :
#            term = lv[1]
#        #end if
        term = lv[0]
        # increment that term's occurrence
        if term in kgSet :
            kgTerms[term] += 1
        else :
            kgTerms[term] = 1
            kgSet.add(term)
        #end if
    #end if


    # stats: collect all the genes, count degrees
    m0 = list()
    m1 = list()
    # check if node matches with a gene regex
    for gt in keepGenes :
        m0.append( re.match(gt, lv[0]) )
        m1.append( re.match(gt, lv[1]) )
    #end loop
    # if a node matched, add to the geneDict
    if any(match != None for match in m0) :
        # add the genes to the dict
        if lv[0] in geneSet :
            geneDict[lv[0]] += 1
        else :
            geneDict[lv[0]] = 1
            geneSet.add(lv[0])
        #end if
    #end if
    if any(match != None for match in m1) :
        # add the genes to the dict
        if lv[1] in geneSet :
            geneDict[lv[1]] += 1
        else :
            geneDict[lv[1]] = 1
            geneSet.add(lv[1])
        #end if
    #end if


    # write to output file
    ofile.write("{1}{0}{2}{0}{3}{0}{4}".format(delim,
        lv[0], lv[1], lv[2], etype))
#end loop
print "There were {:,} lines in the file.".format(count)

# close up shop
ifile.close()
ofile.close()


# Save the list of genes to a file
geneList = list(geneSet)
del geneSet
geneList.sort()
gefile = outname + ".genes.txt"
print "Saving genes to {}".format(gefile)
gef = open(path + gefile, 'wb')
first = True
for g in geneList :
    if first :
        gef.write("{}{}{}".format(g, delim, geneDict[g]))
        first = False
        continue
    #end if
    gef.write("\n{}{}{}".format(g, delim, geneDict[g]))
#end loop
gef.close()

geneDict.clear()
del geneList



######## ######## ######## ########
# Calculate & output some basic graph stats

statfile = outname + '.stat.basic.txt'
print "Writing the stats file: {} ...".format(statfile)

sf = open(path + statfile, 'wb')
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
    sf.write("  {:<20}{:>11,}{:>8.2f}%\n".format( et, eTypes[et], 100 * eTypes[et] / float(numEdges) ))
#end loop
sf.write("\n")

#stats: gene counts
gene_regex = ['ENSG', 'ENSM', 'FBgn', 'GB[45][0123456789]', 'WBGe', 'Y[ABCDEFGHIJKLMNOP][LR][0123456]', 'AT[12345CM]G']
gene_source = ['human', 'mouse', 'fly', 'bee', 'worm', 'yeast', 'grass']




# stats: node misc
nodeList = list(nodeSet)
nodeList.sort()
sf.write("Total Nodes: {:,}\n".format( len(nodeList) ))
sf.write("Node Degree distributions by type:\n")
sf.write("    -  Node Type  -      -# Nodes-     -Mean-  \n")
#of.write("         -  Node Type  -          -# Nodes-     -Mean- -StdDev-  -Min-  -Max- \n")
# first the genes
#of.write(" {:<11} {:>8,} {:9.2f} {:8.2f} {:>5} {:>7}\n".format("genes", len(stats), np.mean(stats), np.std(stats), np.min(stats), np.max(stats) ))
del nodeSet
#del nodeList
nodes.clear()
sf.write("\n")

#-? #numGenes = len(np.unique(Edges[:,1]))
#-? #of.write(" {:<20} {:>11,} {:>11.1f} \n".format("genes", numGenes, numEdges / float(numGenes) ))
#-? # then the rest
#-? for et in typeList :
#-?     numType = len(np.unique(Edges[Edges[:,3] == et][:,0]))
#-?     sf.write(" {:<20} {:>11,} {:>11.1f} \n".format( et, numType, numEdges / float(numType)))
#-? #end loop
#-? sf.write("\n")

# stats: go terms
if "GO_term" in typeList :
    # Extract the go terms
#    goList = nodeList[nodeList[0:3] == "GO:"]
    goList = list(goSet)
    del goSet
    goList.sort()
    numgo = len(goList)

#    print goList

#    # Get the size of each go term
#    goLen = dict()
#    goSet = set()
#    for gt in goList :
#        if gt in goSet :
#            goLen[gt] += 1
#        else :
#            goSet.add(gt)
#            goLen[gt] = 1
#        #end if
#    #end loop

    # Save term lengths to a file
    gofile = outname + ".goList.txt"
    print "Saving GO terms to {}".format(gofile)
    gf = open(path + gofile, 'wb')
#    goList = list(goSet)
#    goList.sort()
    for gt in goList :
        gf.write("{}{}{}\n".format(gt, delim, goTerms[gt]))
    gf.close()

    sf.write("Total GO terms: {:,}\n".format(len(goList)))
    sf.write("Number of GO terms of size ...\n")
    # Count how many terms exist of each length
    lenList = list(goTerms.values())
    lenList.sort()
    for l in range(0,10) :
        a = l * 10
        b = a + 9
#        tempList = lenList[lenList >= a]
#        tempList = tempList[tempList <= b]
        tempList = [1 for x in lenList if (x >= a) & (x <= b)]
        gcount = sum(tempList)
#        gcount = ((lenList >= a) & (lenList <= b)).sum()
        sf.write("  {:>3} -{:>3} :  {:>4}\n".format(a, b, gcount))
    #end loop
    for l in range(2,6) :
        a = l * 50
        b = a + 49
#        tempList = lenList[lenList >= a]
#        tempList = tempList[tempList <= b]
        tempList = [1 for x in lenList if (x >= a) & (x <= b)]
        gcount = sum(tempList)
#        gcount = len(tempList)
#        gcount = ((lenList >= a) & (lenList >= b)).sum()
        sf.write("  {:>3} -{:>3} :  {:>4}\n".format(a, b, gcount))
    #end loop
#    tempList = lenList[lenList >= 300]
    tempList = [1 for x in lenList if (x >= 300)]
    gcount = sum(tempList)
    sf.write("  300+     :  {:>4}\n".format(gcount))
    sf.write('\n')

    del goList
    goTerms.clear()
#end if

# stats: motif terms
if "motif_u5_gc" in typeList :
    mtList = list(mtSet)
    del mtSet
    mtList.sort()
    nummt = len(mtList)

    # Save term lengths to a file
    mtfile = outname + ".mtList.txt"
    print "Saving motif terms to {}".format(mtfile)
    mf = open(path + mtfile, 'wb')
    for mt in mtList :
        mf.write("{}{}{}\n".format(mt, delim, mtTerms[mt]))
    mf.close()

    sf.write("Total motif terms: {:,}\n".format(len(mtList)))
    sf.write("Number of motif terms of size ...\n")
    # Count how many terms exist of each length
    lenList = list(mtTerms.values())
    lenList.sort()
    for l in range(0,10) :
        a = l * 1000
        b = a + 999
        tempList = [1 for x in lenList if (x >= a) & (x <= b)]
        mcount = sum(tempList)
        sf.write("  {:>4}-{:<4} :  {:>4}\n".format(a, b, mcount))
    #end loop
#    for l in range(2,6) :
#        a = l * 50
#        b = a + 49
#        tempList = [1 for x in lenList if (x >= a) & (x <= b)]
#        mcount = sum(tempList)
#        sf.write("  {:>3} -{:>3} :  {:>4}\n".format(a, b, mcount))
#    #end loop
    tempList = [1 for x in lenList if (x >= 300)]
    mcount = sum(tempList)
    sf.write("  10,000+   :  {:>4}\n".format(mcount))
    sf.write('\n')

    del mtList
    mtTerms.clear()
#end if


# stats: Protein Families
if "pfam_domain" in typeList :
    pfList = list(pfSet)
    del pfSet
    pfList.sort()
    numpf = len(pfList)

    # Save term lengths to a file
    pffile = outname + ".pfList.txt"
    print "Saving Protein Families to {}".format(pffile)
    pf = open(path + pffile, 'wb')
    for pt in pfList :
        pf.write("{}{}{}\n".format(pt, delim, pfTerms[pt]))
    pf.close()

    sf.write("Total Protein Family terms: {:,}\n".format(len(pfList)))
    sf.write("Number of Protein Families of size ...\n")
    # Count how many terms exist of each length
    lenList = list(pfTerms.values())
    lenList.sort()
#    for l in range(0,10) :
#        a = l * 10
#        b = a + 9
#        tempList = [1 for x in lenList if (x >= a) & (x <= b)]
#        pfcount = sum(tempList)
#        sf.write("  {:>3} -{:>3} :  {:>4}\n".format(a, b, pfcount))
#    #end loop
#    for l in range(2,7) :
#        a = l * 50
#        b = a + 49
#        tempList = [1 for x in lenList if (x >= a) & (x <= b)]
#        pfcount = sum(tempList)
#        sf.write("  {:>3} -{:>3} :  {:>4}\n".format(a, b, pfcount))
#    #end loop
#    tempList = [1 for x in lenList if (x >= 300)]
#    pfcount = sum(tempList)
#    sf.write("  300+     :  {:>4}\n".format(pfcount))
    for l in range(0,5) :
        a = l * 20
        b = a + 19
        tempList = [1 for x in lenList if (x >= a) & (x <= b)]
        pfcount = sum(tempList)
        sf.write("  {:>3} -{:>3} :  {:>4}\n".format(a, b, pfcount))
    for l in range(1,4) :
        a = l * 100
        b = a + 99
        tempList = [1 for x in lenList if (x >= a) & (x <= b)]
        pfcount = sum(tempList)
        sf.write("  {:>3} -{:>3} :  {:>4}\n".format(a, b, pfcount))
    #end loop
    for l in range(2,5) :
        a = l * 200
        b = a + 199
        tempList = [1 for x in lenList if (x >= a) & (x <= b)]
        pfcount = sum(tempList)
        sf.write("  {:>3} -{:>3} :  {:>4}\n".format(a, b, pfcount))
    #end loop
    tempList = [1 for x in lenList if (x >= 1000)]
    pfcount = sum(tempList)
    sf.write("  1000+    :  {:>4}\n".format(pfcount))
    sf.write('\n')

    del pfList
    pfTerms.clear()
#end if

# stats: KEGG terms
if "KEGG" in typeList :
    kgList = list(kgSet)
    del kgSet
    kgList.sort()
    numkg = len(kgList)

    # Save term lengths to a file
    kgfile = outname + ".kgList.txt"
    print "Saving KEGG terms to {}".format(kgfile)
    kf = open(path + kgfile, 'wb')
    for kg in kgList :
        kf.write("{}{}{}\n".format(kg, delim, kgTerms[kg]))
    kf.close()

    sf.write("Total KEGG terms: {:,}\n".format(len(kgList)))
    sf.write("Number of KEGG terms of size ...\n")
    # Count how many terms exist of each length
    lenList = list(kgTerms.values())
    lenList.sort()
    for l in range(0,10) :
        a = l * 10
        b = a + 9
        tempList = [1 for x in lenList if (x >= a) & (x <= b)]
        kcount = sum(tempList)
        sf.write("  {:>3} -{:>3} :  {:>4}\n".format(a, b, kcount))
    #end loop
    for l in range(2,6) :
        a = l * 50
        b = a + 49
        tempList = [1 for x in lenList if (x >= a) & (x <= b)]
        kcount = sum(tempList)
        sf.write("  {:>3} -{:>3} :  {:>4}\n".format(a, b, kcount))
    #end loop
    tempList = [1 for x in lenList if (x >= 300)]
    kcount = sum(tempList)
    sf.write("  300+     :  {:>4}\n".format(kcount))
    sf.write('\n')

    del kgList
    kgTerms.clear()
#end if


sf.close()




print "\nDone.\n"
####### ####### ####### ####### 