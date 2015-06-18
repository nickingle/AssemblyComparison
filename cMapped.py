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

def plotGraph(d1,d2):
	x = []
	y = []
	
	x = d1.keys()
	y = d1.values()
	
	plt.plot([float(v) for v in x],[float(k) for k in y],'.')
	plt.title('# of Minia Contigs Mapped to Velvet Refs', fontsize = 25)
	plt.ylabel('# of Minia Contigs Mapped to Reference Seq', fontsize = 13)
	plt.yticks(fontsize = 15)
	plt.xlabel('Length of Velvet Sequence', fontsize = 20)
	plt.xticks(fontsize = 15)
	plt.savefig('Contigs_Mapped_SP.png')
	
	x = d2.keys()
	y = d2.values()
	
	plt.plot([float(v) for v in x],[float(k) for k in y], '.')
	plt.title('# of Velvet Contigs Mapped to Minia Refs', fontsize = 20)
	plt.ylabel('# of Velvet Contigs Mapped to Reference Seq', fontsize = 13)
	plt.yticks(fontsize = 15)
	plt.xlabel('Length of Minia Sequence', fontsize = 15)
	plt.xticks(fontsize = 15)
	plt.savefig('Contigs_Mapped_VtoM.png')
	

if __name__ == "__main__":
	input("Press Enter to display results...")
	samfile1 = pysam.AlignmentFile(sys.argv[1], "r")
	samfile2 = pysam.AlignmentFile(sys.argv[2], "r")
	d1 = createDict(samfile1)
	d2 = createDict(samfile2)
	samfile1.close()
	samfile2.close()
	plotGraph(d1,d2)
	
	
	
	
	
				
		
		
		
	
	
	
	
