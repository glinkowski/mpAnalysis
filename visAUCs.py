# ---------------------------------------------------------
# author: Greg Linkowski
# project: Metapath Analysis
#       for the KnowEnG big data center at UIUC
#       funded by the NIH
# 
# Draw area under the curve(s) for predictions
#
# Draw the area under the ROC and the PR curve for a given
#	prediction attempt.
# ---------------------------------------------------------

import matplotlib.pyplot as plt
import visLibrary as vl
import os



####### ####### ####### ####### 
# PARAMETERS

# Paths to network files & sample prediction files
pPath = '../Dropbox/mp/output/pred03-batch-002/'
pFolderPrefix = ''
#pPath = '../Dropbox/mp/output/'
#pFolderPrefix = 'pred01-'

####### ####### ####### ####### 



####### ####### ####### ####### 
# BEGIN MAIN FUNCTION
print ""


# Get the list of subdirectories in pPath
allItems = os.listdir(pPath)
subDirs = [x for x in allItems if os.path.isdir(pPath+x)]


# Create the figure for each valid prediction
for pFolder in subDirs :
	# Skip folders that don't start with desired prefix
	if not pFolder.startswith(pFolderPrefix) :
		continue
	#end if

	if pFolder.endswith('/') :
		pDir = pPath + pFolder
	else :
		pDir = pPath + pFolder + '/'
	#end if




#	# Get the curve statistics from the ranked_genes file
#	FPR, recall, precision, numHid = vl.getAUCstats(pDir)



	# Find all the 'ranked_genes...' files
	fNames = [f for f in os.listdir(pDir) if f.startswith('ranked_genes')]
	fNames.sort()

	for f in fNames :

		# Get the curve statistics from the ranked_genes file
		FPR, recall, precision, numHid = vl.getAUCstats(pDir, f)

#		print FPR
		if FPR[0] == (-1):
			print("  ...skipping...")
			continue
		#end if

		# Calculate (approximate) are under the ROC curve
		areaROC = 0
		for r in recall :
			areaROC += (r / len(recall))
		#end loop

		outName = vl.nameOutputFile(pDir, 'AUC-'+pFolder, 'png')

		# Plot the results
		fig = plt.figure()

		# Plot the ROC curve
		plt.subplot(1, 2, 1)
		plt.plot(FPR, recall)
		plt.plot([0,1], [0,1], 'lightgrey')
		plt.xlabel('False Positive Rate')
		plt.ylabel('True Positive Rate')
		plt.axis([0, 1, 0, 1])

		# Plot the Precision-Recall curve
		plt.subplot(1, 2, 2)
		plt.plot(recall, precision)
		plt.xlabel('Recall')
		plt.ylabel('Precision')
		plt.axis([0, 1, 0, 1])
	#TODO: why is P-R line so jagged?

		# Final touches
		plt.suptitle(pFolder+', concealed = {}'.format(numHid)+
			', ROC area = {:.3}'.format(areaROC))
	#	plt.tight_layout()
		plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=.4, hspace=None)

		# Save the figure
	#	plt.show()
		plt.savefig(pDir+outName)
#		plt.savefig(pDir+'AUC_'+pFolder+'.png')
		plt.close()

	#end loop

#end loop


#TODO: ? Create an overlay of the different methods?


print "\nDone.\n"