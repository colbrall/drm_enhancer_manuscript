# drm_finder.py
# a program for finding DRMs in a FASTA sequence
# takes a FASTA file as input, writes a new file with boxplot statistics
# can count using regular expressions, or by various MOODS settings. For moods,
# specify the background n.t. frequencies you want to use (human, equal, default)
# and the p-value. It can also either output a file with summary stats for each
# DRM, or a BED file with all the enhancers' counts (or both). As long as the
# FASTA file path is first, it defaults to regex and summary file if it can't
# parse the rest of the arguments. P defaults to 0 if it manages to
# parse a MOODS background. It throws an exception if it can't open
# the FASTA file.
#
# Required Arguments:
#	- Path to fasta file
#	- counting method [regex, human, equal, default]
#		if it's a MOODS method:
#		- P value - entered as denominator (e.g. "1024" to indicate 1/1024)
# 	- output file type (summary, all, or both)
#		if all or both:
# 		- path to bed file with same enhancers in same order as FASTA file.
# Example:
#	python drm_finder.py PATH_TO_FASTA method (p_value) output (PATH_TO_BED)
#
# @author Laura Colbran 7/2014; modified 11/2015


import sys
import string
import re
import MOODS
import os
import numpy
import time


# given a sequence string, returns its complement
# only called if writing sequence match file
def complement(line):
        new = ""
        for i in range(len(line)):
            if (line[i] == "A" or line[i] == "a"):
                new += "T"
            elif (line[i] == "T" or line[i] == "t"):
                new += "A"
            elif (line[i] == "G" or line[i] == "g"):
                new += "C"
            else:
                new += "G"
        return new[::-1]


# finds and returns a list of lists of string representations of all DRMs in a sequence.
# can count by regular expressions or by various MOODS methods.
def count(enhlist, c, p_val):
	write_match = False
	index = 0
	seq = ""
	while index <= len(enhlist) - 1 and not enhlist[index].startswith('>'):
		seq += enhlist[index]
		seq = seq.replace('\n','')
		index += 1
	index = 0
	results = []
	pseudocount = 0.001
#	starrmot15
#	caca = [[49,0,100,0,91,0,100,0],[1,83,0,88,0,49,0,94],[49,1,0,12,0,9,0,6],[0,16,0,0,9,41,0,0,]]
	caca = [[0,1,0,1,0,1],[1,0,1,0,1,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]
	me137 = [[0,1,0,1,0,1],[0,0,0,0,0,0],[1,0,1,0,1,0],[0,0,0,0,0,0]]
	gcgc = [[0,0,0,0,0,0],[0,1,0,1,0,1],[1,0,1,0,1,0],[0,0,0,0,0,0]]
	tata = [[0,1,0,1,0,1],[0,0,0,0,0,0],[0,0,0,0,0,0],[1,0,1,0,1,0]]
	matrices = [me137,caca,gcgc,tata]

	if (c == "default"):
		results = MOODS.search(seq,matrices,p_val,both_strands=True)

	elif (c == "human"):
		BG_hum =[0.29508855202553025, 0.20466109233964447, 0.20478482916547036, 0.2954655264693549]
	  	matrices = [MOODS.count_log_odds(matrix,BG_hum,pseudocount) for matrix in matrices]
	  	thresholds = [MOODS.threshold_from_p(matrix, BG_hum, p_val) for matrix in matrices]
	  	results = MOODS.search(seq, matrices, thresholds, convert_log_odds=False, threshold_from_p=False,both_strands=True)

	elif (c == "equal"):
		BG_eq = [0.25,0.25,0.25,0.25]
	  	matrices = [MOODS.count_log_odds(matrix,BG_eq,pseudocount) for matrix in matrices]
	  	thresholds = [MOODS.threshold_from_p(matrix, BG_eq, p_val) for matrix in matrices]
	  	results = MOODS.search(seq, matrices, thresholds, convert_log_odds=False, threshold_from_p=False,both_strands=True)

	else:
		write_match = False
		ga = re.compile(r'GAGA(?:GA)+',re.IGNORECASE)
		ct = re.compile(r'TCTC(?:TC)+',re.IGNORECASE)
		ca = re.compile(r'CACA(?:CA)+',re.IGNORECASE)
		gt = re.compile(r'TGTG(?:TG)+',re.IGNORECASE)
		gc = re.compile(r'GCGC(?:GC)+',re.IGNORECASE)
		cg = re.compile(r'CGCG(?:CG)+',re.IGNORECASE)
		ta = re.compile(r'TATA(?:TA)+',re.IGNORECASE)
		at = re.compile(r'ATAT(?:AT)+',re.IGNORECASE)

		ga_occur = re.findall(ga,seq) + re.findall(ct,seq)
		ca_occur = re.findall(ca,seq) + re.findall(gt,seq)
		gc_occur = re.findall(gc,seq) + re.findall(cg,seq)
		ta_occur = re.findall(ta,seq) + re.findall(at,seq)
		results = [ga_occur, ca_occur, gc_occur, ta_occur]
#	write files containing the matches for each motif as found by MOODS
#   only if boolean set at top of count() is True
	if write_match:
		ga = open('%sga_matches_%d.txt' % (time.strftime("%Y-%m-%d"),1/p_val), 'a')
		ca = open('%sca_matches_%d.txt' % (time.strftime("%Y-%m-%d"),1/p_val), 'a')
		gc = open('%sgc_matches_%d.txt' % (time.strftime("%Y-%m-%d"),1/p_val), 'a')
		ta = open('%sta_matches_%d.txt' % (time.strftime("%Y-%m-%d"),1/p_val), 'a')
		comp = False
		for match in results[0]:
			pos = match[0]
			st = ""
			if (pos < 0):
				comp = True
				pos += len(seq)
			for i in range(len(me137[1])):
				st += seq[pos+i]
			if comp:
				st = complement(st)
				comp = False
			ga.write('%s\n' % st)

		for match in results[1]:
			pos = match[0]
			st = ""
			if (pos < 0):
				comp = True
				pos += len(seq)
			for i in range(len(caca[1])):
				st += seq[pos+i]
			if comp:
				st = complement(st)
				comp = False
			ca.write('%s\n' % st)

		for match in results[2]:
			pos = match[0]
			st = ""
			if (pos < 0):
				comp = True
				pos += len(seq)
			for i in range(len(gcgc[1])):
				st += seq[pos+i]
			if comp:
				st = complement(st)
				comp = False
			gc.write('%s\n' % st)

		for match in results[3]:
			pos = match[0]
			st = ""
			if (pos < 0):
				comp = True
				pos += len(seq)

			for i in range(len(tata[1])):
				st += seq[pos+i]
			if comp:
				st = complement(st)
				comp = False
			ta.write('%s\n' % st)

		ga.close()
		ca.close()
		gc.close()
		ta.close()
	results.append(len(seq))
	return results


# calculates and writes a file containing statistics for the
# given list of drms in each sequence in the fasta file
def stats(list1, list2, file):
	ls = list1
	len_list = list2
	num_enh = len(ls)
	total_drm = 0
	total_len = 0
	u = 0
	l = 0
	total_per_bp = 0
	enh_ls = []
	dense_ls = []
	for i in range(num_enh):
		enh = ls[i]
		leg = len_list[i]
		num_drm = len(enh)
		total_drm += num_drm
		dense_ls.append(float(num_drm)/leg)
		enh_ls.append(num_drm)

	if len(enh_ls) == 0:
		file.write('No DRMs\n')
		file.write('\n~~~~~\n\n')
		return
	ls.sort(key=len)

# 	calculate means and standard deviations
	st_dev_per_enh = numpy.std(enh_ls)
	avg_drm_per_enh = numpy.mean(enh_ls)
	avg_drm_per_bp = numpy.mean(dense_ls)
	st_dev_per_bp = numpy.std(dense_ls)

	file.write('Total Num. DRMs: %d\n' % (total_drm))
	file.write('Mean Num. DRM: %.3f\n' % (avg_drm_per_enh))
	file.write('St. Dev: %.3f\n\n' % (st_dev_per_enh))
#	file.write('Q1: %.3f\t\tMedian: %.3f\t\tQ3: %.3f\n' % (q1_num_drm, med_num_drm, q3_num_drm))
	file.write('Mean Num. DRM/bp: %f\n' % (avg_drm_per_bp))
	file.write('St. Dev: %.3f\n\n' % st_dev_per_bp)
	file.write('\n~~~~~\n\n')


def main():
	try:
		filename = sys.argv[1]
		infile = open(filename, 'r')
		f = str.split(filename,"/")
		f = str.split(f[-1],".")
	except:
		print "Invalid FASTA File Path"
		print "Correct Arguments: PATH_TO_FILE method (p) output (PATH_TO_BED)"
		print "See program comments for details"
		sys.exit()
#	make sure usable arguments are a thing
	method = "regex"
	p_val = 0
	out = "summary"
	bedfile = "NA"
	print "Counting DRMs in %s" % filename
#       out is automatically handled-
#       if it's supposed to be anything other than summary, then it must be
#       at [-2]. If that's not the case, it defaults to summary anyway.
	if (len(sys.argv) > 2):
		method = sys.argv[2]
		out = sys.argv[-2]
		bedfile = sys.argv[-1]
	try:
		p_val = float(1)/int(sys.argv[3])
	except:
		p_val = 0
	if (method != "human" and method != "default" and method != "equal"):
		method = "regex"
		print "Regular Expressions"
	else:
		print "MOODS %s background" % method
		print "p-value: %.6f" % p_val
#       if bedfile is given, it must always be last. If it's not parsible,
#       the program will just output the summary file
	try:
		bed = open(bedfile, 'r')
	except:
		out = "summary"
	if (out != "all" and out != "both"):
		out = "summary"
		print "Output: summary"
	else:
		print "Output: %s" % out
		print "using BEDfile %s" % bedfile
	print "\nCounting DRMs..."
#	count the DRMs
	enhancers = infile.readlines()
	infile.close()
	index = 0
	ga_list = []
	ca_list = []
	gc_list = []
	ta_list = []
	len_list = []
	tot_num = 0
	for line in enhancers:
		if line.startswith('>'):
			count_list = count(enhancers[index+1:], method, p_val)
			ga_list.append(count_list[0])
			ca_list.append(count_list[1])
			gc_list.append(count_list[2])
			ta_list.append(count_list[3])
			len_list.append(count_list[4])
			tot_num += 1
		index += 1

	print "Writing output..."
#	write the complete count file
	if (out == "all" or out == "both"):
		newfile = open('%s_%s_%s.bed' % (time.strftime("%Y-%m-%d"),f[0],method), 'a')
		lines = bed.readlines()
		for i in range(len(ga_list)):
			lines[i] = lines[i].replace('\n', '\t%d\t%d\t%d\t%d\n' % (len(ga_list[i]), len(ca_list[i]), len(gc_list[i]), len(ta_list[i])))
			newfile.write(lines[i])
		newfile.close()
		bed.close()

#	write the summary file
	else:
		newfile = open('%s_%s_%s_summary.txt' % (time.strftime("%Y-%m-%d"),f[0],method), 'a')
		newfile.write('# ')
		for entry in (sys.argv):
			newfile.write('%s ' % entry)
		newfile.write('\n# %s\n'% time.strftime("%c"))
		newfile.write('-------------------------------------\n')
		newfile.write('Num. Enhancers: %d\n\nGAGAGA:\n\n' % (tot_num))
		stats(ga_list, len_list, newfile)
		newfile.write('-------------------------------------\nCACACA:\n\n')
		stats(ca_list, len_list, newfile)
		newfile.write('-------------------------------------\nGCGCGC:\n\n')
		stats(gc_list, len_list, newfile)
		newfile.write('-------------------------------------\nTATATA:\n\n')
		stats(ta_list, len_list, newfile)
		newfile.close()

#	print info about enhancer length
	avg = numpy.mean(len_list)
	min = numpy.amin(len_list)
	max = numpy.amax(len_list)
	print '\nEnhancer Lengths'
	print 'Average: %d' % avg
	print 'Min: %d' % min
	print 'Max: %d\n' % max
	print 'Total Enhancers: %d' % len(len_list)

if __name__ == "__main__":
	main()
