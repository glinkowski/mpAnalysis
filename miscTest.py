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



spath = 'samplesFake/'
sList = ff.getSampleList(spath)
print 'The following samples were found in {}'.format(spath)
print sList
print ''

spath = 'samplesMSIG/'
sList = ff.getSampleList(spath)
print 'The following samples were found in {}'.format(spath)
for item in sList :
    print item
print ''


