import sys
import pysam
import string
import re
import numpy as np
import matplotlib
import matplotlib as plt

def createDict(sfile):
	refLengths = ()
	refLengths = refLengths + sfile.lengths
	
	LengthDict = {}
	
	for readseg in sfile.fetch():
		index = readseg.reference_id
		lengthRef = list(refLengths)[index] # returns the length of the corresponding reference seq that that the read maps to
		if LengthDict.has_key(lengthRef):
			count = list(LengthDict)[lengthRef]
			list(LengthDict)[lengthRef] = count + 1
		else:
			count = 1
			list(LengthDict)[lengthRef] = count
	
	return LengthDict

def plotGraph(d):
	x = []
	y = []
	
	x = d.keys()
	y = d.values()
	
	plt.plot(x,y,'.')
	plt.title('# of Minia Contigs Mapped to Velvet Refs', fontsize = 25)
	plt.ylabel('# of Minia Contigs Mapped to Reference Seq', fontsize = 20)
	plt.yticks(fontsize = 15)
	plt.xlabel('Length of Velvet Sequence', fontsize = 20)
	plt.xticks(fontsize = 15)
	plt.savefig('Contigs_Mapped_SP.png')

if __name__ == "__main__":
	input("Press Enter to display results...")
	samfile = samfile = pysam.AlignmentFile(sys.argv[1], "r")
	d = createDict(samfile)
	samfile.close()
	plotGraph(d)
	
	
	
	
	
				
		
		
		
	
	
	
	
