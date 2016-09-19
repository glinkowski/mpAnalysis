Repo: mpAnalysis
# Metapath Analysis of Genomic Data
## author: Greg Linkowski
for the KnowEnG big data center at UIUC, 
funded by the NIH

##Hypothesis:
We have access to a multitude of data regarding how genes are related to each other. Yet, for specific cancers and other diseases, we still don't know the specific mechanism that connects the affected genes. Here, I use an analysis of meta-paths -- specific patterns of edge types in a graph -- to approximate those connections. This allows us to describe the connections between the genes in question in terms of weighted combinations of the data we have. Further work will allow us to drill down deeper into that data to uncover biologically relevant details describing those connections.

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
2. Extract alternative features by running preFeatureTerms.py and/or preFeatureNeighbor.py
3. Depending on data source, may need to pre-format the gene sets using preFormatSamples.py or preFormatAchilles.py

######Part 02: characterize the gene sets
1. Run mpChar01.py to extract set-similarity features and create cross-validation folds
2. Run mpChar02.py to rank features and genes by similarity to the set
3. (future) Run mpChar03.py to provide analysis across folds and sets


##Description of the scripts in this repository

###Part 01: pre-process the network  
Must be done before any further work. Allows fast statistics generation on network.
* preProcessing -- script which applies the steps below to the network
* preProcCont -- script to continue pre-processing if interrupted
* preProcFuncs -- library of functions used during pre-processing
* preFeatureNeighbor -- creates feature vectors describing node degree rather than direct connections
* preFeatureTerms -- creates feature vectors where membership to individual terms is the only consideration
* preFormatSamples -- extracts the sets from the jumbled data output from the KnowEnG network
* preFormatAchilles -- fixes file names to match the expected input

#####Steps performed during pre-processing:

1. fix known typos, exclude specified edge & node types, normalize edge weights and apply threshold if desired
2. remove non-gene nodes, drawing new edges directly between genes
3. save the modified network edge file, list of genes, and node-index dictionary
4. define nodes as high/med/low -degree, for use in random sample selection
5. create the primary path matrices for each edge type in network
6. create and save the metapath matrices from the primaries


###Part 02: (previous version) select important features for a set and rank remaining genes by similarity
* mpPredict04a -- for a batch of samples, create the PathSim-based similarity vectors
* mpPredict04a2 -- for a batch of samples, create the z-score-based similarity vectors
* mpPredict04g -- select relevant features and predict gene membership ranking using clustering and selected classifier
* mpPredict04i -- select features and rank genes using iterative/voting with Lasso approach
* mpPredict04x -- compile prediction data for a batch of samples, find AUC scores
* mpLibrary -- functions used as part of the metapath prediction process, will replace mpFindFuncs

###Part 02 (current version) Characterize gene sets through analysis of meta-paths
* mpChar01 -- create the z-score similarity vectors, and create cross-validation folds
* mpChar02 -- rank genes by similarity and detail important / useful features
* (future) mpChar03 -- analysis of results, AUC scores of folds


