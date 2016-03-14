Repo: 15-Falls
# Metapath Analysis of Genomic Data

##How to use the contained scripts
####Current Approach:
######Part 01: pre-process the network
0. Given an edge file, create a .keep.txt and optional .correct.txt file. Examples can be seen in the ./network folder.
  1. The keep file indicates which edge and gene species to keep, as well as indicating which edge types are "indirect." That is, if an edge begins or ends at a node that is not a gene. For edge types, the name of the edge type is given. For kept genes, the Python regular expression for the first four characters of the gene is given.
  2. The correct file contains known mis-spellings on the left, with the correction on the right.
1. Run preProcessing.py
  1. NOTE: It is assumed the files are named: (name).edge.txt, (name).keep.txt, (name).correct.txt
  2. Set the parameters at the beginning of the file
    1. ename, epath: the name of the network (ie: (name) from above) and the path to the files
    2. mpDepth: the maximum metapath length to calculate, up to maximum of 4 (NOTE: time grows exponentially)
  3. Run the script: preProcessing.py
  4. Continue an interrupted build with preProcCont.py
    1. If the network processing is interrupted, or for example, if you previously used mpDepth=2 and now want to calculate the 3-step metapaths, set the calcStep parameter to the level you wish to continue calculating.


##Description of the scripts in this repository
###Part 00: modify the network 
* mkSubNtwk00 -- select certain edge and node types to create a sub-network, these functions have been subsumed by the preProcessing script
* mkToy00 -- trim the network, keeping only a certain percent of the nodes (create a smaller network for testing)

###Part 01: pre-process the network  
Must be done before any further work. Allows fast statistics generation on network.
* preProcessing -- script which applies the steps below to the network
* preProcCont -- script to continue pre-processing if interrupted
* preProcFuncs -- library of functions used during pre-processing

#####Steps performed during pre-processing:

1. fix known typos, exclude specified edge & node types, normalize edge weights and apply threshold if desired
2. remove non-gene nodes, drawing new edges directly between genes
3. save the modified network edge file, list of genes, and node-index dictionary
4. define nodes as high/med/low -degree, for use in random sample selection
5. create the primary path matrices for each edge type in network
6. create and save the metapath matrices from the primaries


###Part 02: collect metapaths for a given sample
* mpFindPaths01 -- output metapath statistics for a single sample
* mpFindPaths02 -- same as above, but use node binning to select random samples with degree distributions similar to the test sample
* mpFindFuncs -- functions used by the mpFindPaths scripts
* sampleMatrix -- functions used for the process of node-binning
* mpPredict01 -- for a single sample, partition and apply a classifier based on metapath analysis in an attempt to predict hidden genes
* mpLibrary -- functions used as part of the metapath prediction process, will replace mpFindFuncs


**Previous Approach**: collect path info using DFS
* mpDFSWalk00 -- use DFS on edge file, tracking paths between nodes

