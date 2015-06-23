import sys
import pysam
import string
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.ioff()

def createDict(sfile):
	refLengths = ()
	refLengths = refLengths + sfile.lengths
	
	LengthDict = {}
	
	for readseg in sfile.fetch():
		index = readseg.reference_id
		lengthRef = refLengths[index] # returns the length of the corresponding reference seq that that the read maps to
		if lengthRef in LengthDict:
			count = LengthDict[lengthRef]
			LengthDict[lengthRef] = count + 1
		else:
			count = 1
			LengthDict[lengthRef] = count
	
	return LengthDict
	
def getOutliers(sfile):
	refID = ()
	refID = refID + sfile.references
	refLengths = ()
	refLengths = refLengths + sfile.lengths
	
	LengthDict = {}
	
	for readseg in sfile.fetch():
		index = readseg.reference_id
		lengthRef = refID[index] # returns the length of the corresponding reference seq that that the read maps to
		if lengthRef in LengthDict:
			count = LengthDict[lengthRef]
			LengthDict[lengthRef] = count + 1
			if LengthDict[lengthRef] > 10:
				print('Reference ID: {0}'.format(lengthRef))
		else:
			count = 1
			LengthDict[lengthRef] = count
	

def plotGraphs(mList, names, n):
	x = []
	y = []
	axlist = list(names)
	num = int(n/2)
	fig, axes = plt.subplots(nrows=num, ncols=2)
	
	for k in range(0, num):
		d = list(mList)[k] # get dictionary from list to map
		x = np.array(list(d.keys()))
		y = np.array(list(d.values()))
		
		axes[k,0].plot(x, y, '.')
		gName = list(names)[k]
		title = '{0} : Minia Contigs mapped to Velvet'.format(gName)
		axes[k,0].set_title(title, fontsize=10)
		axes[k,0].set_ylabel('# of Minia Contigs Mapped to Reference Seq', fontsize=7)
		axes[k,0].set_xlabel('Length of Velvet Sequence', fontsize=10)
		
		d = list(mList)[k+1]
		x = np.array(list(d.keys()))
		y = np.array(list(d.values()))
		axes[k,1].plot(x, y, '.')
		title = '{0} : Velvet Contigs mapped to Minia Sequences'.format(gName)
		axes[k,1].set_title(title, fontsize=10)
		axes[k,1].set_ylabel('# of Velvet Contigs Mapped to Reference Seq', fontsize=7)
		axes[k,1].set_xlabel('Length of Minia Sequence', fontsize=10)
	
	fig.subplots_adjust(top=2)
	plt.tight_layout()
	fig.savefig('Contigs_Mapped.png')
	print('Saved figure to Contigs_Mapped.png')
	plt.close(fig)
	

if __name__ == "__main__":	
	###### From main.py for taking in multiple input files
	n = int(sys.argv[1])
	cm_array = []
	inputnum = 3*n
	numOfFiles = n*2
	counter = 0
	samfilelist = list()
	cmList = list()
	names = list()
	
	for x in range(2,(inputnum+2)):
		counter += 1
		if counter == 3:
			names.append(str(sys.argv[x]))
			counter = 0
		else:
			samfile = pysam.AlignmentFile(str(sys.argv[x]), "r")
			samfilelist.append(samfile) # [x-2-numOfnames] is location of 
			d = createDict(samfile)
			getOutliers(samfile)
			cmList.append(d)
			samfile.close()
	numOfFilesI = int(numOfFiles)		
	plotGraphs(cmList, names, numOfFiles)
	

	
				
		
		
		
	
	
	
	
