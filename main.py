import sys
import pysam
import string
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Method for creating dictionaries: 
# MisMatch Dictionary: mm_d = [qname, MD Tag] for each line of SAM
# CigarString Dictionary: cs_d = [qname], CigarString] for each line of SAM
def createDicts(sfile):
	mm_d = {}
	cs_d = {}
	for readseg in sfile.fetch():
		qname = readseg.qname
		try:
			md_string = readseg.get_tag("MD")
			
		except KeyError:
			md_string = "None"
			continue
		mm_d[qname] = md_string
		
		ct = readseg.cigartuples
		cs_d[qname] = ct
	
	return mm_d, cs_d

# Method for breaking down MD tag by mismatches:
# seperates(split) the string by the numbers of matched ntides (A decimal num) 
# from the mismatched ntides (chars A,C,G,T) 
def inter_md(md):
	md_list = re.split("([A-Z]+)", md)
	return md_list

# Method that returns the Match Percent of the MD tag by taking in a
# MD list returned by inter_md method and counting the number of matches
# and comparing it to Matches + MisMatches
def match_percent(l):
	count = 0
	matches = 0
	
	for t in l:
		if t.isdigit():
			num = int(t)
			matches += num
			count += num
		else:
			count += 1
	#print matches
	#print count
	percent = float(matches)/float(count)
	
	return percent
	
# Method that calculates the Percent Identity for a given cigar string
# Also counts the number of other signs in the CIGAR, but only returns 
# a number that is a percent of matches and it ignores the clipped ends
def calc_pi(cig):
	matches = 0
	insertions = 0
	deletions = 0
	skips = 0
	softclip = 0
	hardclip = 0
	padding = 0
	seqmatch = 0
	seqmismatch = 0
	
	for l in cig:
		if l[0] == 0:
			matches += l[1]
		elif l[0] == 1:
			insertions += l[1]
		elif l[0] == 2:
			deletions += l[1]
		elif l[0] == 3:
			skips += l[1]
		elif l[0] == 4:
			softclip += l[1]
		elif l[0] == 5:
			hardclip += l[1]
		elif l[0] == 6:
			padding += l[1]
		elif l[0] == 7:
			seqmatch += l[1]
		elif l[0] == 8:			
			seqmismatch += l[1]
		else:
			continue
	
	pid = float(matches)/(float(matches)+float(skips)+float(insertions)+float(softclip)+float(hardclip))
	
	return pid
	
# Method that calculates the Percent of the sequence that is matched and covered
def calc_perMatch(cig, md):
	# Break down Cigar String ########################### 
	matches = 0
	insertions = 0
	deletions = 0
	skips = 0
	softclip = 0
	hardclip = 0
	padding = 0
	seqmatch = 0
	seqmismatch = 0
	
	for l in cig:
		if l[0] == 0:
			matches += l[1]
		elif l[0] == 1:
			insertions += l[1]
		elif l[0] == 2:
			deletions += l[1]
		elif l[0] == 3:
			skips += l[1]
		elif l[0] == 4:
			softclip += l[1]
		elif l[0] == 5:
			hardclip += l[1]
		elif l[0] == 6:
			padding += l[1]
		elif l[0] == 7:
			seqmatch += l[1]
		elif l[0] == 8:			
			seqmismatch += l[1]
		else:
			continue
	###############################################################
	
	# Interpret MD tag ############################################
	count = 0
	md_matches = 0
	
	if md == 0:
		percent = float(matches)/(float(matches)+float(skips)+float(insertions)+float(softclip)+float(hardclip))
		return percent
	else:	
		for t in md:
			if t.isdigit():
				num = int(t)
				md_matches += num
				count += num
			else:
				count += 1
	#####################################################################
	
	percent = float(md_matches)/(float(matches)+float(skips)+float(insertions)+float(softclip)+float(hardclip))
	percent = percent*100
	return percent

# Method for creating dictionary that includes the qname and percent match 
# using both the cigar string and the md tag
def createPMDic(csd, mdd):
	pm_d = {}
	
	for s in csd:
		name = list(csd)[num]
		cig = csd[name]
		md_string = mdd[name]
		if md_string != "None":
			mdl = inter_md(md_string)
		else:
			mdl = 0
		pm = calc_perMatch(cig,mdl)
		pm_d[name] = pm
	
	return pm_d	
	
def createPMarray(csd, mdd):
	pm_a = []
	num = 0
	
	for s in csd:
		name = list(csd)[num]
		cig = csd[name]
		md_string = mdd[name]
		if md_string != "None":
			mdl = inter_md(md_string)
		else:
			mdl = 0
		pm = calc_perMatch(cig,mdl)
		pm_a.append(pm)
		num += 1
	
	return pm_a

def printHighMatches(pma):
	count = 0
	
	for m in pma:
		if m >= 98:
			count += 1
		else:
			continue
	print "# of Alignments w/ a Match Percentage over 98 percent = ", count			

def graphPM(pma):
	printHighMatches(pma)
	plt.hist(pma, range = (min(pma), 99))
	plt.ylabel('# of Alignments')
	plt.xlabel('Percent Match')
	plt.title('Count of Percent Matches')
	plt.show()
					

# Method that creates a series for NumPy to interpret and plot using the 
# qname and its MD match percent
def createMDSeries(mdd):
	md_sd = {}
	num = 0
	n_count = 0
	
	for tag in mdd:
		name = list(mdd)[num]
		md_string = mdd[name]
		if md_string != "None":
			mdl = inter_md(md_string)
			md_match = match_percent(mdl) 
			md_sd[num] = md_match
			num += 1
		else:
			n_count += 1
			num += 1
	
	md_series = pd.Series(md_sd)
	return md_series, n_count			

# Method that creates a series for NumPy to interpret and plot using the
# qname and its cigar string percent identity
def createCigSeries(csd):
	cs_sd = {}
	num = 0
	
	for t in csd:
		name = list(csd)[num]
		cs = csd[name]
		pi = calc_pi(cs)
		num += 1
		cs_sd[num] = pi
		
	
	cs_series = pd.Series(cs_sd)
	return cs_series

									
# MAIN METHOD ##################################################################
# ##############################################################################
if __name__ == "__main__":
	samfile = pysam.AlignmentFile(sys.argv[1], "r")
	md_dict, cig_dict = createDicts(samfile)
	pm_array = createPMarray(cig_dict, md_dict)
	graphPM(pm_array)
	
	
	
	samfile.close()
	
# Finished createPercent Match array to build histogram graph to show coverage
# percentages
	
'''
cseries = createCigSeries(cig_dict)
mdseries, num = createMDSeries(md_dict)

mdseries.plot(title = "MD Tag Percent Match")
cseries.plot(title = "Cigar String Percent Match")
plt.show()
'''	
	
	
	
	
	
	
	
	
