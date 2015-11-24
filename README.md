Repo: 15-Falls
# Metapath Analysis of Genomic Data

* mpfuncs -- frequently-used functions

**Current Approach**: recreate the network as edge matricies, 
		 building a master metapath matrix, so that testing a
		 sample takes minutes, not days.

**Part 01**: pre-process the network  
Must be done before any further work. Allows fast statistics generation on network.
* preProc00 -- fix known typos, collect basic network stats
* preProc01 -- normalize edge weights per edge type
         node-index dictionary used for fast edge lookups give node name
* preProc02 -- replace indirect edges with direct gene-gene edges
         create path matrix for each edge type (gene-by-gene)
* preProc03 -- final pre-processing step: calculate the metapaths
         up to a certain length, save in an output file


**Part 02**: collect metapaths for a given sample
* demo01 -- proof-of-concept, outputs file with comparison
		 to random samples of similar size


**Previous Approach**: collect path info using DFS
* mpDFSWalk00 -- use DFS on edge file, tracking paths between nodes

