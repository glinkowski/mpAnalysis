# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# Pre-Processing of the network
#
# This is to continue creating metapaths at higher levels
#   after the preProcessing.py script has been run.
# Outline:
# ---------------------------------------------------------

import preProcFuncs as pp
import time



####### ####### ####### ####### 
# PARAMETERS

# The network to use and directory path
ename = 'all_v3beta_g2e9t0'
epath = '../Dropbox/mp/networks/'
#ename = 'fakeNtwk00_g2e3t10'
#epath = 'networks/'

# The metapath level to calculate
calcStep = 3

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

tstart = time.time()

# Load the primary matrices to memory
print("\nLoading the primary path matrices")
primList = pp.readPrimaryMatrices(epath, ename)
#primNames = list(primList.keys())
#primNames.sort()
#print primNames
#print primList[0]

mpPath = epath + ename + "_MetaPaths/"
# Calculate the requested metapaths
print("Creating metapaths of length {}".format(calcStep))
if calcStep == 1 :
    pp.createMPLengthOne(primList, mpPath)
if calcStep == 2 :
    pp.createMPLengthTwo(primList, mpPath)
elif calcStep == 3 :
    pp.createMPLengthThreeFast(primList, mpPath)
elif calcStep == 4 :
    pp.createMPLengthFour(primList, mpPath)
#end if


print("\nDone.\n")