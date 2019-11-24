# @description: Merge bam files of multiple sequencing runs of same library prep 

#### SAMPLES IN SEQ RUN NEED TO HAVE THE SAME SAMPLE NAME 

import re
import os
import sys
import getopt

#### IMPORT FILELIST NEED TO CONTAIN SAMPLE NAME, INPUT PATH FIRST BAMFILE, INPUT PATH SECOND BAMFILE

### ADD PATH TO SAMBAMBA HERE:
sambamba_path = "XXX"

try:
	opts, args = getopt.getopt(sys.argv[1:], 'f:o:', ['filelist', 'output'])
except getopt.GetoptError:
	usage()
	sys.exit(2)

for opt, arg in opts:
	if opt in ('-f', '--filelist'):
		import_filelist = arg
	elif opt in ('-o', '--output'):
		output_path = arg
	else:
		usage()
		sys.exit(2)

with open(import_filelist) as f1:
	for line in f1:
		line = line.rstrip()
		columns = line.split("\t")
		sample = columns[0]
		input_path_run1 = columns[1]
		input_path_run2 = columns[2]
		
		sample_run1 = input_path_run1 + "/" + sample + "/mapping/" + sample + "_sort.bam" 
		sample_run2 = input_path_run2 + "/" + sample + "/mapping/" + sample + "_sort.bam"
		
		os.system("cd " + output_path + "; " + "mkdir " + sample)
		os.system("cd " + output_path + "/" + sample + "; " + "mkdir mapping")
		output_merge = output_path + "/" + sample + "/mapping/" + sample + "_merge.bam"
		output_sort = output_path + "/" + sample + "/mapping/" + sample + "_merge.sorted.bam"
		
		print(output_merge) 
		
		os.system("java -jar $PICARD MergeSamFiles I=" + sample_run1 + " I=" + sample_run2 + " O=" + output_merge)
		os.system(sambamba_path + " sort " + output_merge)
		os.system(sambamba_path + " index " + output_sort)
		#os.system("rm " + output_merge)