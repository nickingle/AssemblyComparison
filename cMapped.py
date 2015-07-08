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
	refLengths = sfile.lengths
	refNames = sfile.references
	refNameL = []
	sumOfLength = 0
	contigCount  = 0
	tuplelength = len(refLengths)
	#print('Length of length tuple: {0}'.format(str(tuplelength)))
	
	#for l in refLengths:
		#sumOfLength = sumOfLength + int(l)  
	#print('Sum of the Reference Lengths: {0}'.format(str(sumOfLength)))
	
	LengthDict = {}
	
	for readseg in sfile.fetch():
		contigCount = contigCount + 1
		index = readseg.reference_id
		refname = refNames[index]
		#if refname in refNameL:
		#	continue
		#else:
		#	refNameL.append(refname)
		#	print('refName: {0}'.format(str(refname)))
		lengthRef = refLengths[index] # returns the length of the corresponding reference seq that that the read maps to
		if lengthRef in LengthDict:
			count = LengthDict[lengthRef]
			LengthDict[lengthRef] = count + 1
			if count > 5000:
				print("Reference Sequence: {0}".format(refname))
		else:
			count = 1
			LengthDict[lengthRef] = count
			
	#print('Contig Count: {0}'.format(str(contigCount)))
	return LengthDict

def plotGraphs(mList, names):
	x = []
	y = []
	axlist = list(names)
	fig, (ax0, ax1) = plt.subplots(ncols=2)
	
	d = list(mList)[0] # get dictionary from list to map
	x = list(d.keys())
	y = list(d.values())
		
	ax0.plot(x, y, '.')
	gName = list(names)[0]
	title = '{0} : Minia Contigs mapped to Velvet'.format(gName)
	ax0.set_title(title, fontsize=10)
	ax0.set_ylabel('# of Minia Contigs Mapped to Reference Seq', fontsize=7)
	ax0.set_xlabel('Length of Velvet Sequence', fontsize=10)
		
	x = []
	y = []
	d = list(mList)[1]
	x = list(d.keys())
	y = list(d.values())
	ax1.plot(x, y, '.')
	title = '{0} : Velvet Contigs mapped to Minia Sequences'.format(gName)
	ax1.set_title(title, fontsize=10)
	ax1.set_ylabel('# of Velvet Contigs Mapped to Reference Seq', fontsize=7)
	ax1.set_xlabel('Length of Minia Sequence', fontsize=10)
		
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
			#print('Name: {0}'.format(str(sys.argv[x])))
		else:
			samfile = pysam.AlignmentFile(str(sys.argv[x]), "r")
			samfilelist.append(samfile) # [x-2-numOfnames] is location of 
			d = createDict(samfile)
			cmList.append(d)
			samfile.close()
	numOfFilesI = int(numOfFiles)
	size = len(cmList)
	print("{0} items in the list".format(size))		
	plotGraphs(cmList, names)


	
				
		
		
		
	
	
	
	
