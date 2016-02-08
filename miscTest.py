# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# This is a test script to call and verify different functions
#
# ---------------------------------------------------------

import mpFindFuncs as ff
import preProcFuncs as pp
import numpy as np


#spath = 'samplesFake/'
#sList = ff.getSampleList(spath)
#print 'The following samples were found in {}'.format(spath)
#print sList
#print ''
#
#spath = 'samplesMSIG/'
#sList = ff.getSampleList(spath)
#print 'The following samples were found in {}'.format(spath)
#for item in sList :
#    print item
#print ''


ename = 'fakeNtwk01'
epath = 'networks/'
kfile = ename + '.keep.txt'
print "\nReading in the network:", ename
print "    reading the keep file", kfile
hGenes, keepGenes, loseGenes, keepEdges, indirEdges, thresh = pp.readKeepFile(epath+kfile)
print 'Keep genes: {}'.format(keepGenes)
print 'Ignore genes: {}'.format(loseGenes)
print 'Keep edges: {}'.format(keepEdges)
print 'Indirect: {}'.format(indirEdges)
print 'Threshold = {}'.format(thresh)
print 'Humans: {}'.format(hGenes)



## Identify unique edge names in a network
#print ''
#fe = open('networks/all_v1.edge.txt', 'rb')
#edgeSet = set()
#for line in fe :
#	line = line.rstrip()
#	lv = line.split('\t')
#
#	if lv[3] not in edgeSet :
#		edgeSet.add(lv[3])
##end loop
#edgeList = list(edgeSet)
#edgeList.sort()
#print 'Edges in all_v1: {}'.format(edgeList)



print ''
ename = 'all_v1'
epath = 'networks/'
kfile = ename + '.keep.txt'
# Read in the keep file
#	if not there, ask to have it created
print "\nReading in the network:", ename
print "    reading the keep file", kfile
geneHuman, keepGenes, loseGenes, keepEdges, indirEdges, thresh = pp.readKeepFile(epath+kfile)
print 'Keep genes: {}'.format(keepGenes)
print 'Ignore genes: {}'.format(loseGenes)
print 'Keep edges: {}'.format(keepEdges)
print 'Indirect: {}'.format(indirEdges)
print 'Threshold = {}'.format(thresh)
print 'Humans: {}'.format(hGenes)
