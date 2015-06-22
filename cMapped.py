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

def plotGraphs(mList, names, n):
	x = []
	y = []
	axlist = list(names)
	num = int(n/2)
	fig, axes = plt.subplots(nrows=num, ncols=2)
	
	for k in range(0, num-1):
		d = list(mList)[k] # get dictionary from list to map
		x = np.array(list(d.keys()))
		y = np.array(list(d.values()))
		
		axes[k,0].plot(x, y, '.')
		gName = list(names)[k]
		title = '{0} : Minia Contigs mapped to Velvet'.format(gName)
		axes[k,0].set_title(title, fontsize=14)
		axes[k,0].set_xlabel('# of Minia Contigs Mapped to Reference Seq', fontsize=10)
		axes[k,0].set_ylabel('Length of Velvet Sequence', fontsize=10)
		
		d = list(mList)[k+1]
		x = np.array(list(d.keys()))
		y = np.array(list(d.values()))
		axes[k,1].plot(x, y, '.')
		title = '{0} : Minia Contigs mapped to Velvet'.format(gName)
		axes[k,1].set_title(title, fontsize=14)
		axes[k,1].set_xlabel('# of Minia Contigs Mapped to Reference Seq', fontsize=10)
		axes[k,1].set_ylabel('Length of Velvet Sequence', fontsize=10)
	
	fig.suptitle("Number of Contigs Mapped to Reference Sequences by Length", fontsize=18)
	fig.savefig('Contigs_Mapped.png')
	plt.close(fig)
	
	
'''	
	x = d1.keys()
	y = d1.values()
	
	plt.plot([float(v) for v in x],[float(k) for k in y],'.')
	plt.title('# of Minia Contigs Mapped to Velvet Refs', fontsize = 20)
	plt.ylabel('# of Minia Contigs Mapped to Reference Seq', fontsize = 13)
	plt.yticks(fontsize = 15)
	plt.xlabel('Length of Velvet Sequence', fontsize = 15)
	plt.xticks(fontsize = 15)
	plt.savefig('Contigs_Mapped_SP.png')
	plt.close("all")
	
	x = d2.keys()
	y = d2.values()
	
	plt.plot([float(v) for v in x],[float(k) for k in y], '.')
	plt.title('# of Velvet Contigs Mapped to Minia Refs', fontsize = 20)
	plt.ylabel('# of Velvet Contigs Mapped to Reference Seq', fontsize = 13)
	plt.yticks(fontsize = 15)
	plt.xlabel('Length of Minia Sequence', fontsize = 15)
	plt.xticks(fontsize = 15)
	plt.savefig('Contigs_Mapped_VtoM.png')
'''	

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
	
	for x in range(2,(inputnum+1)):
		counter += 1
		if counter == 3:
			names.append(str(sys.argv[x]))
			counter = 0
		else:
			samfile = pysam.AlignmentFile(str(sys.argv[x]), "r")
			samfilelist.append(samfile) # [x-2-numOfnames] is location of 
			d = createDict(samfile)
			cmList.append(d)
			samfile.close()
	numOfFilesI = int(numOfFiles)		
	plotGraphs(cmList, names, numOfFiles)
	
'''	
	input("Press Enter to display results...")
	samfile1 = pysam.AlignmentFile(sys.argv[1], "r")
	samfile2 = pysam.AlignmentFile(sys.argv[2], "r")
	d1 = createDict(samfile1)
	d2 = createDict(samfile2)
	samfile1.close()
	samfile2.close()
	plotGraph(d1,d2)	
'''	
	
				
		
		
		
	
	
	
	
