# ----------------------------------------------------
# EXPERIMENT 16-1
# Goal: Create a specific sub-network from given network
# ----------------------------------------------------

import gfuncs2 as gf
import numpy as np
import re



######## ######## ######## ########
# PARAMETERS

edgepath = "edgefiles/"
edgefile = "all.edge.txt"
#edgefile = "hsa_dghmw.edge.txt"
delim = "\t"

outpath = "edgefiles/"
outfile = "h_g_t0.edge.txt"
# filename: species_edges_threshold.edge.txt

# threshold: keep only weights above X (range [0,1])
threshold = 0

# species               gene designation
include_human = True   # ENSG
include_mouse = False   # ENSM
include_fly = False     # FBgn
include_bee = False     # GB[45][0123456789]
include_worm = False    # WBGe
include_yeast = False   # Y[ABCDEFGHIJKLMNOP][LR][0123456]
include_grass = False   # AT[12345CM]G

# edge types                designation
include_GOTerms = True     # GO_term
include_ProtDomain = False  # pfam_domain
include_motif = False       # motif_u5_gc
include_homology = False    # prot_homol
include_KEGG = False        # KEGG
include_biogrid = False     # PPI_BioGRID
include_dip = False         # PPI_DIP
include_mint = False        # PPI_MINT
include_intact = False      # PPI_IntAct
include_SCoexp = False      # STRING_coexpression
include_SCoOccur = False    # STRING_cooccurrence
include_SDatabase = False   # STRING_database
include_SExp = False        # STRING_experimental
include_SFusion = False     # STRING_fusion
include_SNeihb = False      # STRING_neighborhood
include_STextM = False      # STRING_textmining

######## ######## ######## ########



######## ######## ######## ########
# IMPLEMENT PARAMETERS

edge_types = ['GO_term', 'pfam_domain', 'motif_u5_gc', 'prot_homol', 'KEGG', 'PPI_BioGRID',  'PPI_DIP', 'PPI_MINT', 'PPI_IntAct', 'STRING_coexpression', 'STRING_cooccurrence', 'STRING_database', 'STRING_experimental', 'STRING_fusion', 'STRING_neighborhood', 'STRING_textmining']
edge_select = [include_GOTerms, include_ProtDomain, include_motif, include_homology, include_KEGG, include_biogrid, include_dip, include_mint, include_intact, include_SCoexp, include_SCoOccur, include_SDatabase, include_SExp, include_SFusion, include_SNeihb, include_STextM]

#gene_regex = ['ENSG', 'GO:0', 'AT[12345CM]G']
gene_regex = ['ENSG', 'ENSM', 'FBgn', 'GB[45][0123456789]', 'WBGe', 'Y[ABCDEFGHIJKLMNOP][LR][0123456]', 'AT[12345CM]G']
gene_select = [include_human, include_mouse, include_fly, include_bee, include_worm, include_yeast, include_grass]

#keepEdges = np.where(edge_select, edge_types)
#keepEdges = [edge_types if edge_select]

keepEdges = list()
for i in range(0, len(edge_types)) :
    if edge_select[i] :
        keepEdges.append(edge_types[i])
#end loop

#print keepEdges

keepGenes = list()
for i in range(0, len(gene_regex)) :
    if gene_select[i] :
        keepGenes.append(gene_regex[i])
#end loop

#print keepGenes

######## ######## ######## ########



######## ######## ######## ########
# BEGIN MAIN FUNCTION

print ""

print "keeping edges of type:", keepEdges
print "keeping genes of type:", keepGenes

print "Reading in the graph ..."
Edges, Nodes = gf.readGraphData(edgepath + edgefile, delim)

newEdges = np.empty( [0,4], dtype='a25')

#print np.shape(newEdges)
#print np.shape(Edges)

for et in keepEdges :
    newEdges = np.vstack( (newEdges, Edges[Edges[:,3] == et]) ) 
#end loop
del Edges


print "Example: "
print newEdges[0:9,:]


#print newEdges


keepRow = list()
#for i in range(0,1) :
for i in range(0, np.shape(newEdges)[0]) :

    # check columns 0 & 1 for genes
#    isGene0 = isGene1 = False
    m0 = list()
    m1 = list()
    for gt in gene_regex :
        m0.append( re.match(gt, newEdges[i,0]) )
        m1.append( re.match(gt, newEdges[i,1]) )
    #end loop
#    print m0
    isGene0 = any(match != None for match in m0)
    isGene1 = any(match != None for match in m1)
#    print isGene0

    # if both columns are genes
    if (isGene0 and isGene1) :
#        print "both genes"

        m0 = list()
        m1 = list()
        for gt in keepGenes :
            m0.append( re.match(gt, newEdges[i,0]) )
            m1.append( re.match(gt, newEdges[i,1]) )
        #end loop
        doKeep0 = any(match != None for match in m0)
        doKeep1 = any(match != None for match in m1)

        # keep if both genes are in keepGenes
        if (doKeep0 and doKeep1) :
            keepRow.append(True)
        else :
            keepRow.append(False)
        #end if
    # if col 0 is a gene, but not col 1
    elif (isGene0 == True) :
#        print "gene 0"

        m0 = list()
        for gt in keepGenes :
            m0.append( re.match(gt, newEdges[i,0]) )
        #end loop
        doKeep0 = any(match != None for match in m0)
        # keep if gene is in keepGenes
        if (doKeep0 == True) :
            keepRow.append(True)
        else :
            keepRow.append(False)
        #end if
    # if col 1 is a gene, but not col 2
    elif (isGene1 == True) :
#        print "gene 1"

        m1 = list()
        for gt in keepGenes :
            m1.append( re.match(gt, newEdges[i,1]) )
        #end loop
        doKeep1 = any(match != None for match in m1)
        # keep if gene is in keepGenes
        if (doKeep1 == True) :
            keepRow.append(True)
        else :
            keepRow.append(False)
        #end if
    # if neither col are genes
    else :
        keepRow.append(False)
    #end if
#end loop

#print keepRow

#print newEdges[0:9,:]
#print keepRow[0:9]

# Doing this the hard way ...
keepIndex = list()
for k in range(0, len(keepRow)) :
    if keepRow[k] :
        keepIndex.append(k)
#end loop
finalEdges = newEdges[keepIndex]
print np.shape(newEdges), np.shape(finalEdges)
del newEdges

#print keepIndex[0:9]
#print finalEdges[0:9,:]


#print np.shape(newEdges)
#print np.shape(finalEdges)


print "Saving new sub-network to file: {}".format(outfile)
#np.savetxt(outpath + outfile, finalEdges)
fout = open(outpath + outfile, "wb")
for r in range(0, np.shape(finalEdges)[0]) :
    fout.write("{0}{4}{1}{4}{2}{4}{3}\n".format(finalEdges[r,0], finalEdges[r,1], finalEdges[r,2], finalEdges[r,3], delim))
#end loop
fout.close()




#TODO: normalize edge weights and apply the threshold

print "\nDone.\n"