# ----------------------------------------------------
# author: Greg Linkowski
# project: metapaths for KnowEnG
# 
# Metapths, step 01
#
# Normalize weights per edge type.
# ----------------------------------------------------




####### ####### ####### ####### 
# PARAMETERS

ename = 'all-v1'
path = '../../edgefiles'

infile = ename + '.edge.txt'
outfile = ename + '.edge_norm.txt'
#outfile = ename + '.gene-only.txt'
delim = '\t'


#indirect=['GO_term', 'motif_u5_gc', 'pfam_domain'] 
#direct = ['prot_homol']

####### ####### ####### ####### 




####### ####### ####### ####### 
# BEGIN MAIN FUNCTION