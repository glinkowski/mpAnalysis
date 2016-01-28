# -*- coding: utf-8 -*-
"""
# ---------------------------------------------------------
# author: Shengliang Dai
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
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

"""


import mpFindFuncs as ff
import pandas as pd
import numpy as np
import random


def sampleMatrix(sname):

    spath = '../samples/'
    opath = '../sampleMatrix/'
    saperatedGenes = '../saperated genes/'
    high = ff.readSampleFiles(saperatedGenes + 'High1nodeBining', False, False)
    med = ff.readSampleFiles(saperatedGenes + 'Med1nodeBining', False, False)
    low = ff.readSampleFiles(saperatedGenes + 'Low1nodeBining', False, False)

    # read a gene list and convert it into numbers
    sampleDict = ff.readFileAsIndexDict('../networks/geneList.txt')
    
   


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

    # do another 99 simulation and make a matrix 
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
    SampleMatrix.to_csv(opath + sname + '_sampleMatrix.txt', sep = '\t', index = False, header = False)
    
    
sampleNames = ff.getSampleList('../samples/')
# read a dictory, select out file names
sampleList = ff.getSampleList('../samples/')

for samples in sampleList:
    sampleMatrix(samples)