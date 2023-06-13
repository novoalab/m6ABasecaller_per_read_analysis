#!/usr/bin/env python 
__author__ = 'luca.cozzuto@crg.eu'

#MODULES
import sys
import re
import optparse
import collections
import pprint
import copy 
import subprocess
import textwrap 
from collections import Counter
from collections import defaultdict


#BODY FUNTIONS
def options_arg():
	usage = "usage: %prog -i <input bam file> -o <OUTPUT PREFIX>"
	parser = optparse.OptionParser(usage=usage)
	parser.add_option('-i', '--input', help='Input bam file', dest="input" )
	parser.add_option('-q', '--min_qual', default = 0, type= int, help='Mininum quality', dest="minqual" )
	parser.add_option('-o', '--output', help='Ouput prefix File', dest="wotus" )
	parser.add_option('-c', '--chromosome', default = "all", help='Chromosome or all', dest="chr" )
	parser.add_option('-s', '--simple', default = "no", help='yes or no', dest="simple" )
	(opts,args) = parser.parse_args()
	if opts.input:pass
	else: parser.print_help()
	return (opts)
	
def __main__ ():
	parsefile(opts.input, opts.wotus, opts.minqual, opts.chr, opts.simple)


#AUXILIAR MODULES
def parsefile(file, oprefix, minqual, chr, simple):
	import pysam

	fw = open(oprefix, "w")

	outText = "readname N_mod ref_name is_reverse start end readlen alnlen is_secondary mod_list\n"
	pp = pprint.PrettyPrinter(indent=4)
	# prepare out files
	samfile = pysam.AlignmentFile(file, "rb")
	if (chr == "all"):
		chr = None
	for read in samfile.fetch(chr):
		if not read.is_unmapped:
			readname = read.query_name
			readlen = read.query_length
			alnlen = read.query_alignment_length
			pos = read.get_reference_positions()
			start = min(pos)
			end = max(pos)

			#Check strand:
			strand = read.is_reverse

			#end = read.query_name
			secondary = read.is_secondary
			#sequence = read.query_sequence 
			quality = read.query_qualities
			#ref_id = read.reference_id
			ref_name = samfile.get_reference_name(read.reference_id)
			#print(ref_name)
			start = read.reference_start
			# Read the cigar content
			#cigar_string = read.cigarstring
			if (quality is not None and sum(quality) >= minqual):
				aln = dict(read.get_aligned_pairs())
				muts = []

				for pos_read,qual in enumerate(quality):
					# got a modification
					if (qual >= minqual):
						muts = addPos(muts, aln, pos_read, simple)
						#mut_pos =  getPos(aln, pos_read)
						#muts.append(str(mut_pos))
				
				mod_n = len(muts)
				if (mod_n>0):
					lastf = ",".join(muts)
				else:
					lastf = "NA"
				outText = outText + readname + " " + str(mod_n) + " " + ref_name + " "  + str(strand) + " " + str(start) + " " + str(end) + " " + str(readlen) + " " + str(alnlen) + " " + str(secondary) + " " + lastf +"\n"
				#else:
				#	outText = outText + readname + " 0"  + " " + ref_name +  + " " + str(readlen) + " NA" +"\n"				
			else:
				lastf = "NA"
				mod_n = 0
				outText = outText + readname + " " + str(mod_n) + " " + ref_name + " " + str(strand) + " " + str(start) + " " + str(end) + " " + str(readlen) + " " + str(alnlen) + " " + str(secondary) + " " + lastf +"\n"

	fw.write(outText)
	fw.close()
							
	return

def addPos(muts, dictionary, index, simple):
	pos = str(dictionary[index])
	indexA = index
	indexB = index
	num = 0
	offset = 1
	if (simple == "no"):
		if (pos == "None"):
			posA = "None"
			posB = "None"
			while(True):
				num = num + 1
				indexA = index + num
				indexB = index - num
				if indexA in dictionary:
					posA = str(dictionary[indexA])
				if indexB in dictionary:
					posB = str(dictionary[indexB])
				if (posA != "None" or posB != "None"):
					if (posA != "None"):
						pos = str(int(posA)+offset) + "(+" + str(num) + ")"
					elif (posB != "None"):
						pos = str(int(posB)+offset) + "(-" + str(num) + ")"
				break
		else:
			pos = str(int(pos) + offset)
		muts.append(str(pos))
	else: 
		if (pos != "None"):
			pos = str(int(pos) + offset)
			muts.append(str(pos))
		
	return muts
		


#Calling
opts = options_arg()
if opts.input and opts.wotus:
	__main__()

