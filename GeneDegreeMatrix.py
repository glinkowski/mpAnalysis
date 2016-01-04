# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 09:59:11 2015

@author: shengliangdai
"""
import sys
import pandas as pd

def degree():
    outname = 'toy_hsa_c'
    path = '../networks/'
    infile = 'toy2_hsa.edge.txt'
	
    outfile = outname + '.degree.txt'

    if (path+outfile) == (path+infile) :
	    print "ERROR: Input and output are same file: {}".format(path+outfile)
	    sys.exit()
    types = ['GO_term', 'motif_u5_gc', 'pfam_domain', 'prot_homol']
    edges = pd.read_table(path + infile, '\t', header = None)
    
    edges = edges.rename(columns = {0:'nodes', 1:'genes', 2:'weight', 3:'type'})

    dictionary = {}
    dictionaryGO = {}
    dictionaryMotif = {}
    dictionaryPfam = {}
    dictionaryProt = {}
    
    for i in range(len(edges.nodes)):
        if edges.nodes[i][0:4] == "ENSG":
            if edges.nodes[i] not in dictionary:
                dictionary[edges.nodes[i]] = 1
            else:
                dictionary[edges.nodes[i]] += 1
            if edges.type[i] == types[0]:
                if edges.nodes[i] not in dictionaryGO:
                    dictionaryGO[edges.nodes[i]] = 1
                else:
                    dictionaryGO[edges.nodes[i]] += 1
            if edges.type[i] == types[1]:
                if edges.nodes[i] not in dictionaryMotif:
                    dictionaryMotif[edges.nodes[i]] = 1
                else:
                    dictionaryMotif[edges.nodes[i]] += 1
            if edges.type[i] == types[2]:
                if edges.nodes[i] not in dictionaryPfam:
                    dictionaryPfam[edges.nodes[i]] = 1
                else:
                    dictionaryPfam[edges.nodes[i]] += 1
            if edges.type[i] == types[3]:
                if edges.nodes[i] not in dictionaryProt:
                    dictionaryProt[edges.nodes[i]] = 1
                else:
                    dictionaryProt[edges.nodes[i]] += 1
    
    for i in range(len(edges.genes)):
        if edges.genes[i][0:4] == "ENSG":
            if edges.genes[i] not in dictionary:
                dictionary[edges.genes[i]] = 1
            else:
                dictionary[edges.genes[i]] += 1   
            if edges.type[i] == types[0]:
                if edges.genes[i] not in dictionaryGO:
                    dictionaryGO[edges.genes[i]] = 1
                else:
                    dictionaryGO[edges.genes[i]] += 1
            if edges.type[i] == types[1]:
                if edges.genes[i] not in dictionaryMotif:
                    dictionaryMotif[edges.genes[i]] = 1
                else:
                    dictionaryMotif[edges.genes[i]] += 1
            if edges.type[i] == types[2]:
                if edges.genes[i] not in dictionaryPfam:
                    dictionaryPfam[edges.genes[i]] = 1
                else:
                    dictionaryPfam[edges.genes[i]] += 1
            if edges.type[i] == types[3]:
                if edges.genes[i] not in dictionaryProt:
                    dictionaryProt[edges.genes[i]] = 1
                else:
                    dictionaryProt[edges.genes[i]] += 1
 
    degreeMatrix = pd.DataFrame.from_dict(dictionary.items())
    degreeGO = pd.DataFrame.from_dict(dictionaryGO.items())
    degreeMotif = pd.DataFrame.from_dict(dictionaryMotif.items())
    degreePfam = pd.DataFrame.from_dict(dictionaryPfam.items())
    degreeProt = pd.DataFrame.from_dict(dictionaryProt.items())
    #degreeMatrix.to_csv(path + outfile, sep = '\t', index = False, header = False)
    #degreeGO.to_csv(path + 'go.txt', sep = '\t', index = False, header = False)
    #degreeMotif.to_csv(path + 'motif.txt', sep = '\t', index = False, header = False)
    #degreePfam.to_csv(path + 'pfam.txt', sep = '\t', index = False, header = False)
    #degreeProt.to_csv(path + 'Prot.txt', sep = '\t', index = False, header = False)
    matrix = pd.DataFrame.merge(degreeMatrix, degreeGO, how = 'left', on = 0)
    matrix = pd.DataFrame.merge(matrix, degreeMotif, how = 'left', on = 0)
    matrix = pd.DataFrame.merge(matrix, degreePfam, how = 'left', on = 0)
    matrix = pd.DataFrame.merge(matrix, degreeProt, how = 'left', on = 0)
    matrix = matrix.rename(columns={list(matrix)[0]:'genes', list(matrix)[1]:'All degree', list(matrix)[2]:'GO_term', list(matrix)[3]:'motif_u5_gc', list(matrix)[4]:'pfam_domain', list(matrix)[5]:'prot_homol'})
    matrix.to_csv(path + 'degreeMatrix.txt', sep = '\t', index = False, header = False)
    
degree()    
    
    
    
        
