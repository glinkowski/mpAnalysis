# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 19:17:44 2016

@author: shengliangdai
"""

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

import sys
import numpy as np
import pandas as pd

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
        High[0].to_csv(path + 'High' + str(columnNumber) + outname, sep = '\t', index = False, header = False) 
        Low[0].to_csv(path + 'Low' + str(columnNumber) + outname, sep = '\t', index = False, header = False)
        Med[0].to_csv(path + 'Med' + str(columnNumber) + outname, sep = '\t', index = False, header = False)
    #High = degreeMatrix[[0, 1]].sort(1, ascending = False).head(high) 
    #Low = degreeMatrix[[0, 1]].sort(1, ascending = False).tail(low)
    #Med = degreeMatrix[[0, 1]].sort(1, ascending = False).reset_index().loc[high:(len(degreeMatrix)-low-1), :]
    #print High[0], Low[0], Med[0]
nodeBining(0.1, 0.3, 'degreeMatrix.txt')
