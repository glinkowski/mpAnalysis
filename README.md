Repo: 15-Falls
# Metapath Analysis of Genomic Data

**Current Approach**: recreate the network as edge matricies, 
		 building a master metapath matrix, so that testing a
		 sample takes minutes, not days. (Or hours, depending
		 on the underlying hardware.)

##Part 00: modify the network 
* mkSubNtwk00 -- select certain edge and node types to create a sub-network, these functions have been subsumed by the preProcessing script
* mkToy00 -- trim the network, keeping only a certain percent of the nodes (create a smaller network for testing)

##Part 01: pre-process the network  
Must be done before any further work. Allows fast statistics generation on network.
* preProcessing -- script which applies the steps below to the network
* preProcCont -- script to continue pre-processing if interrupted
* preProcFuncs -- library of functions used during pre-processing

####Steps performed during pre-processing:

1. fix known typos, exclude specified edge & node types, normalize edge weights and apply threshold if desired
2. remove non-gene nodes, drawing new edges directly between genes
3. save the modified network edge file, list of genes, and node-index dictionary
4. define nodes as high/med/low -degree, for use in random sample selection
5. create the primary path matrices for each edge type in network
6. create and save the metapath matrices from the primaries


##Part 02: collect metapaths for a given sample
* mpFindPaths01 -- output metapath statistics for a single sample
* mpFindPaths02 -- same as above, but use node binning to select random samples with degree distributions similar to the test sample
* mpFindFuncs -- functions used by the mpFindPaths scripts
* sampleMatrix -- functions used for the process of node-binning
* mpPredict01 -- for a single sample, partition and apply a classifier based on metapath analysis in an attempt to predict hidden genes
* mpLibrary -- functions used as part of the metapath prediction process, will replace mpFindFuncs


**Previous Approach**: collect path info using DFS
* mpDFSWalk00 -- use DFS on edge file, tracking paths between nodes

