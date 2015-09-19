# ----------------------------------------------------
# author: Greg Linkowski
# project: metapaths for KnowEnG
# 
# Metapths, step 01
#
# Convert the network to only gene-gene edges.
# ----------------------------------------------------




####### ####### ####### ####### 
# PARAMETERS

ename = 'all-v1'
path = '../../edgefiles'

infile = ename + 'edge.txt'
outfile = ename + '.gene-only.txt'
delim = '\t'


indirect=['GO_term', 'motif_u5_gc', 'pfam_domain'] 
direct = ['prot_homol']

####### ####### ####### ####### 




####### ####### ####### ####### 
# BEGIN MAIN FUNCTION