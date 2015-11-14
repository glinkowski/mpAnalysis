Repo: 15-Falls
# Metapath Analysis of Genomic Data

mpfuncs -- frequently-used functions

Part 01: pre-processing the network
       Must be done before any further work. Allows fast statistics generation on network.
 mp00 -- fix known typos, collect basic network stats
 mp01 -- normalize edge weights per edge type
         node-index dictionary used for fast edge lookups give node name
 mp02 -- replace indirect edges with direct gene-gene edges
         create path matrix for each edge type (gene-by-gene)
 mp03 -- final pre-processing step: calculate the metapaths
         up to a certain length, save in an output file

