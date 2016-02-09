# -*- coding: utf-8 -*-
"""

# ---------------------------------------------------------
# author: ShengliangDai
# project: Metapath Analysis
#		for the KnowEnG big data center at UIUC
#		funded by the NIH
# 
# functions used to Find Paths in the processed network
#
# These functions were created to aid in the reading of a
#	network, creation of random samples based on high/med/low, and sample matrices
#	
#   
#
# Functions provided:
#	degree()
#	nodeBining(thresholdHigh, thresholdLow, infile)
#	sampleMatrix(sname)
#	writeNodeBinFiles(ntwkPath, ntwkName, geneList, edgeArray, geneHead, binThresholds)
#	degreeMatrix(edgeArray, geneHead)
#	

# ---------------------------------------------------------

"""

'''
# Creat degree matrices for original networks
#
#  Load a network and compare the genes in samples with the genes in original network. 
#   Count the degree for each gene.
#   Build a matrix to store the degree and the degree for each type of gene.
#   
#   
# Outline:
#   1) Read in a network, for example, toy2_hsa.
#   2) Classify genes from the network
#   3) Count how many genes in the given sample are high, median or low degree genes.
#   4) Straitified sampling genes with the same distribution as 3) from original network.
#   5) Map gene names into indices to save space.
#   6) Output the results to a file
# ---------------------------------------------------------
'''


import sys
import pandas as pd
import mpFindFuncs as ff
import numpy as np
import random
import os.path

'''
README
writeNodeBinFiles:
The goal for this function is to separeate genes by their degree to high/med/low bins.
Input:
    the edge array of the original network
    the gene head such as ENSG
    the gene list of the original network

Output:
    write high/med/low bins into .txt files for future sample generationg
    
'''

def writeNodeBinFiles(ntwkPath, ntwkName, geneList, edgeArray, geneHead, binThresholds):
    
    
    outputFolder = ntwkPath + ntwkName + '/'
    
    # If folder doesn't exist, create it
    if not os.path.exists(outputFolder):
        os.makedirs(outputFolder)
    
    
    degreeMatrix2, headers = degreeMatrix(edgeArray, geneHead)
    
    degreeMatrix2 = pd.DataFrame(degreeMatrix2)
    thresholdHigh = binThresholds[0]
    thresholdLow = 1 - binThresholds[1]
    high = int(len(degreeMatrix2)*thresholdHigh)
    
    low = int(len(degreeMatrix2)*thresholdLow) 
    
    
    # Extract each column of types of the matrix, in order to make code flexible to every matrix
    for columnNumber in range(1, len(list(degreeMatrix2))):
        High = degreeMatrix2[[0, columnNumber]].sort(columnNumber, ascending = False).head(high) 
        Low = degreeMatrix2[[0, columnNumber]].sort(columnNumber, ascending = False).tail(low)
        Med = degreeMatrix2[[0, columnNumber]].sort(columnNumber, ascending = False).reset_index().loc[high:(len(degreeMatrix2)-low-1), :]
        High[0].to_csv(outputFolder + 'High' + str(columnNumber)  + '.txt', sep = '\t', index = False, header = False) 
        Low[0].to_csv(outputFolder + 'Low' + str(columnNumber) + '.txt', sep = '\t', index = False, header = False)
        Med[0].to_csv(outputFolder + 'Med' + str(columnNumber) + '.txt', sep = '\t', index = False, header = False)
    
    return 
    
    
'''
README
degreeMatrix:
The goal for this function is to creat a degree matrix from an edge array.
Input:
    the edge array of the original network
    the gene head such as ENSG

Output:
    an np ndarray matrix of degree matrix
    
'''

def degreeMatrix(edgeArray, geneHead):
    edges = pd.DataFrame(edgeArray)
    edges = edges.rename(columns = {0:'nodes', 1:'genes', 2:'weight', 3:'type'})
    types = np.unique(edges.type) #['typeA' 'typeB' 'typeC']
    print types
    all_degree = {}
    # Creat dictionaries based on the number of types, count the degree of genes based on gene types.
    dictName = {}
    for j in range(len(types)):
        dictionaryName = 'dictionary_' + str(types[j])
        dictName[dictionaryName] = {}
    print dictName
    for i in range(len(edges.nodes)):
        if edges.nodes[i][0:4] in geneHead:
            # For the first column of edge matrix, count all degree
            if edges.nodes[i] not in all_degree:
                all_degree[edges.nodes[i]] = 1
            else:
                all_degree[edges.nodes[i]] += 1
            
            # For the first column of edge matrix, count the degree based on different gene types
            if edges.nodes[i] not in dictName['dictionary_' + str(edges.type[i])]:
                dictName['dictionary_' + str(edges.type[i])][edges.nodes[i]] = 1
            else:
                dictName['dictionary_' + str(edges.type[i])][edges.nodes[i]] += 1
                
        if edges.genes[i][0:4] in geneHead:    
            # For the second column of edge matrix, count all degree
            if edges.genes[i] not in all_degree:
                all_degree[edges.genes[i]] = 1
            else:
                all_degree[edges.genes[i]] += 1
                
            # For the second column of edge matrix, count the degree based on different gene types 
            if edges.genes[i] not in dictName['dictionary_' + str(edges.type[i])]:
                dictName['dictionary_' + str(edges.type[i])][edges.genes[i]] = 1
            else:
                dictName['dictionary_' + str(edges.type[i])][edges.genes[i]] += 1
                          
    degreeMatrix = pd.DataFrame.from_dict(all_degree.items())
    for j in range(len(types)):
        tempDegreeMatrix = pd.DataFrame.from_dict(dictName['dictionary_' + str(types[j])].items())
        print tempDegreeMatrix
        degreeMatrix = pd.DataFrame.merge(degreeMatrix, tempDegreeMatrix, how = 'left', on = 0)
    #Fill missing values
    degreeMatrix = degreeMatrix.fillna(0)
    
    headers = ['Genes', 'All_degree',] + list(types)
    
    return degreeMatrix.values, headers

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
    
'''
README
Node Binning:
The goal for this function is to separate out high and low-degree nodes.
Input:
    the list of genes in the network
    the matrix created with the function degreeMatrix.py
Output:
    three lists (or sets) of genes: high, med, low
    
To run the code 
Set threshold for high, low degrees (thresholdHigh, threshholdLow)
Matrix should have the format [Genes, All degree, GO_term, ...]
'''


def nodeBining(thresholdHigh, thresholdLow, infile):
    outname = 'nodeBining.txt'
    path = '../networks/'
    
    
    degreeMatrix = pd.read_table(path + infile, '\t', header = None)
    #Fill missing values
    degreeMatrix = degreeMatrix.fillna(0)
    
    high = int(len(degreeMatrix)*thresholdHigh)
    low = int(len(degreeMatrix)*thresholdLow) 
    
    
    # Extract each column of types of the matrix, in order to make code flexible to every matrix
    for columnNumber in range(1, len(list(degreeMatrix))):
        High = degreeMatrix[[0, columnNumber]].sort(columnNumber, ascending = False).head(high) 
        Low = degreeMatrix[[0, columnNumber]].sort(columnNumber, ascending = False).tail(low)
        Med = degreeMatrix[[0, columnNumber]].sort(columnNumber, ascending = False).reset_index().loc[high:(len(degreeMatrix)-low-1), :]
        #High[0].to_csv(path + 'High' + str(columnNumber) + outname, sep = '\t', index = False, header = False) 
        #Low[0].to_csv(path + 'Low' + str(columnNumber) + outname, sep = '\t', index = False, header = False)
        #Med[0].to_csv(path + 'Med' + str(columnNumber) + outname, sep = '\t', index = False, header = False)
    #High = degreeMatrix[[0, 1]].sort(1, ascending = False).head(high) 
    #Low = degreeMatrix[[0, 1]].sort(1, ascending = False).tail(low)
    #Med = degreeMatrix[[0, 1]].sort(1, ascending = False).reset_index().loc[high:(len(degreeMatrix)-low-1), :]
    #print High[0], Low[0], Med[0]
    return High[0].values

'''
# Creat sample index matrices for specific samples
#
#  Load all sample files and compare the genes in samples with the genes in original network. 
#   Classify genes in original network into three groups based on their degrees: 
#   high degree genes, median degree genes, low degree genes.
#   
#   
# Outline:
#	1) Read in the sample file
#	2) Classify genes from the network
#	3) Count how many genes in the given sample are high, median or low degree genes.
#	4) Straitified sampling genes with the same distribution as 3) from original network.
#   5) Map gene names into indices to save space.
#   6) Output the results to a file
# ---------------------------------------------------------
'''

def sampleMatrix(ntwkPath, ntwkName, geneList, spath, sname):

    #spath = '../samples/'
    #spath = 'samplesFake/'
    #opath = '../sampleMatrix/'
    saperatedGenes = ntwkPath + ntwkName + '/' 
    high = ff.readSampleFiles(saperatedGenes + 'High1', False, False)
    med = ff.readSampleFiles(saperatedGenes + 'Med1', False, False)
    low = ff.readSampleFiles(saperatedGenes + 'Low1', False, False)

    # read a gene list and convert it into numbers
    sampleDict = ff.readFileAsIndexDict(geneList)
    
   


    print "Finding metapaths in sample: {}".format(sname)
    sampGenes = ff.readSampleFiles(spath + sname, True, True)
    sampGenes = pd.Series(sampGenes)
    High = pd.Series(high)
    Med = pd.Series(med)
    Low = pd.Series(low)

    high = sampGenes.isin(High)
    med = sampGenes.isin(Med)
    low = sampGenes.isin(Low)

    print "This sample contains {} high degree nodes".format(np.sum(high))
    print "This sample contains {} med degree nodes".format(np.sum(med))
    print "This sample contains {} low degree nodes".format(np.sum(low))

    #Randomly select genes from network in order to match the distribution of the sample
    rows = random.sample(High.index, np.sum(high))
    HighSample = High[rows]
    HighSample = pd.Series(ff.convertToIndices(HighSample, sampleDict))
    
    rows = random.sample(Med.index, np.sum(med))
    MedSample = Med[rows]
    MedSample = pd.Series(ff.convertToIndices(MedSample, sampleDict))

    rows = random.sample(low.index, np.sum(low))
    LowSample = Low[rows]
    LowSample = pd.Series(ff.convertToIndices(LowSample, sampleDict))
    


    randomSample = HighSample.append(MedSample)
    randomSample = randomSample.append(LowSample)
    randomSample = randomSample.reset_index(drop=True)
    SampleMatrix = randomSample

    for i in range(99):
        i += 1

        rows = random.sample(High.index, np.sum(high))
        HighSample = High[rows]
        HighSample = pd.Series(ff.convertToIndices(HighSample, sampleDict))
        
        rows = random.sample(Med.index, np.sum(med))
        MedSample = Med[rows]
        MedSample = pd.Series(ff.convertToIndices(MedSample, sampleDict))
        
        rows = random.sample(low.index, np.sum(low))
        LowSample = Low[rows]
        LowSample = pd.Series(ff.convertToIndices(LowSample, sampleDict))
        

        randomSample = HighSample.append(MedSample)
        randomSample = randomSample.append(LowSample)
        randomSample = randomSample.reset_index(drop=True)
        SampleMatrix = pd.concat((SampleMatrix,randomSample), axis = 1)

    #print SampleMatrix
    SampleMatrix = SampleMatrix.T
    #SampleMatrix.to_csv(opath + sname + '_sampleMatrix.txt', sep = '\t', index = False, header = False)
    return SampleMatrix.values







        
    
    
