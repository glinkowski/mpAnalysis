# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# Find all the metapaths present in all available samples
#
# This allows for comparison of metapaths across a set of
#   samples. It reads in all the smples in a folder, then
#   outputs the z-Scores as a spreadsheet.
# Outline:
#
    # get a list of all samples in folder
    # create array to store z-Scores
    # create the metapath dict
    # SLOW METHOD ??
    #   for each sample
    #       read the file(s)
    #       create the random samples
    #       calculate statistics
    #       add to the z-Scores array
    # output to a file
# ---------------------------------------------------------

import mpFindFuncs as ff
import time

## REMOVE: The following is for verification.
#import random
#random.seed(42)



####### ####### ####### ####### 
# PARAMETERS

# The network to use and directory path
#ename = 'fakeNtwk01_g2e3t10'
ename = 'fakeNtwk01_g3e4t1'
epath = 'networks/'

# The path to the sample files
spath = 'samplesFake/'

# Where to store the output
oname = 'find02-' + ename
opath = 'outputFake/'

# How many random samples to examine
numRand = 100

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION
print ""

tstart = time.time()


# 1) Get the list of metapaths for this network

print "Checking what paths are available ..."
pathDict = ff.readKeyFile(epath, ename)
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)


# 2) Get the list of all samples in the folder

print "Checking what samples are in {}".format(spath)
sampList = ff.getSampleList(spath)
#print sampList
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)


# 3) Collect the statistics

# Calculate the z-Scores
print "Calculating the statistics"
zScores = ff.createZScoreMatrix(epath, ename,
    pathDict, spath, sampList, numRand)
#print zScores
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)


# 4) Output the data to a file
print "Saving data to ..."
outFile = ff.writeOutputMultSamples(opath, oname, ename,
    sampList, pathDict, zScores)
print "    ... {}".format(outFile)
print "    --elapsed time: {:.3} (s)".format(time.time()-tstart)


#TODO: output the left-out genes?
#   This could be done as a separate script
#   with less overhead


print "\nDone.\n"