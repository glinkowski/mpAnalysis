# ----------------------------------------------------
# author: Greg Linkowski
# project: metapaths for KnowEnG
# 
# Metapths, step 02
#
# Create a network consisting only of gene-gene edges.
#	Replace indirect edges with direct edges. Along
#	the way, create new edges reflective of the
#	strength of membership for that term.
# For example, 3-4 new GO_term edge types to represent
#	higher import of sharing a smaller term.
# ----------------------------------------------------




####### ####### ####### ####### 
# PARAMETERS

ename = 'all-v1'
path = '../../edgefiles'

infile = ename + '.edge.txt'
outfile = ename + '.edge_norm.txt'
#outfile = ename + '.gene-only.txt'
delim = '\t'


indirect=(['GO_term', 'motif_u5_gc', 'pfam_domain', 'KEGG', 
	'PPI_BioGRID'])
direct = (['prot_homol', 'STRING_experimental', 'STRING_coexpression',
	'STRING_textmining', 'STRING_neighborhood', 'PPI_IntAct', 
	'STRING_database', 'STRING_cooccurrence', 'STRING_fusion', 
	'PPI_MINT', 'PPI_DIP'])

# for example:
go_cutoffs = [30, 100, 300]

####### ####### ####### ####### 




####### ####### ####### ####### 
# BEGIN MAIN FUNCTION