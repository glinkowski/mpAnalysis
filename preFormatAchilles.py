


import mpLibrary as mp
import os.path
from os import listdir
import sys
import numpy as np





####### ####### ####### ####### 
# PARAMETERS

sPath = '../Dropbox/mp/samplesAchilles1st/subset02/'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION

mp.verifyDirectory(sPath, False, False)

print('Renaming the samples in {}'.format(sPath))

# Get list of all text files in folder
fNames = [f for f in listdir(sPath) if f.endswith('.txt')]
# if verbose :
# 	print (fNames)

for origName in fNames :

#	# strip the '.txt'
#	origName = origName[:-4]

	# strip the gene count at end (and '.txt')
	modName = origName.split('-')[0]

	# move the '_DN' or '_UP' to end
	ov = modName.split('_')
	newName = ov[0]
	for i in range(2, len(ov)) :
		newName = newName + '_' + ov[i]
	newName = newName + '_' + ov[1] + '.txt'

	# copy the file to the new name
	os.rename(sPath + origName, sPath + newName)
#end loop


print('\nDone.\n')