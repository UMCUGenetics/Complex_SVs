#!/usr/bin/python
# this script is used to run the manta-flag script over multiple Manta outputs

import sys
import getopt
import re
import os
import shutil

try:
    opts, args = getopt.getopt(sys.argv[1:], 'f:h:o:d:l:p:', ['filelist', 'help', 'output', 'distance', 'limit','parents'])
except getopt.GetoptError:
    usage()
    sys.exit(2)


for opt, arg in opts:
    if opt in ('-h', '--help'):
        usage()
        sys.exit(2)
    elif opt in ('-f', '--filelist'):
        filelist = arg
    elif opt in ('-o', '--output'):
        output_folder = arg
    elif opt in ('-d', '--distance'):
        dist = arg
    elif opt in ('-l', '--limit'):
        limit = arg
    elif opt in ('-p', '--parents'):
        parents = arg
    else:
        usage()
        sys.exit(2)

### REPLACE XXX WITH INPUT DIRECTORY
project_dir = "XXX"
processed_dir = "XXX"

filter_script = "python " + input_dir + "/Flag_Manta.py"

temp_dir = output_folder + "temp/"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

if not os.path.exists(temp_dir):
    os.makedirs(temp_dir)

def run_flag_script( mode, input_file, filter, distance, label, output_file, threshold = limit, filter_script = filter_script ):
	if mode != "AC":
		command = filter_script + " -m " + mode + " -i " + input_file + " -f " + filter + " -d " + distance + " -l " + label + " > " + output_file
	
	else:
		command = filter_script + " -i " + input_file + " -t " + threshold + " -l " + label + " > " + output_file
	
	print command
	os.system(command)
	return 
	
with open(filelist) as f1:
	for line in f1:
		line = line.rstrip()
		print line
		
		if not line.startswith("#"):
			columns = line.split("\t")
			sample = columns[0]			
			name = columns[1]
			folder = columns[2]
			
			print "### Start SV flagging for " + sample
			
			input_manta =  processed_dir + folder + "/structuralVariants/manta/" + sample + "/results/variants/diploidSV.vcf.gz"
			
			delly_file = project_dir + "SV/Data/Delly/" + sample + "_delly.vcf"
				
			AC_script = "python " + project_dir + "Scripts/SV/Manta_filter_AltCov.py"			
								
			VCF_1000G = project_dir + "Common_data/1000G/ALL.wgs.integrated_sv_map_v2.20130502.svs.vcf"
			
			VCF_DGV = project_dir + "Common_data/DGV/GRCh37_hg19_variants_2016-05-15_filtered.txt"
			
			manta_DB_filtered = processed_dir + folder + "/structuralVariants/manta/" + sample + "/results/variants/diploidSV.dbfilter.vcf.gz"
					
			exclusion_file = project_dir + "Common_data/Lumpy/ceph18.b37.lumpy.exclude.2014-01-15.bed"		
						
			out1 = temp_dir + name + "_manta_out1.vcf"
			out2 = temp_dir + name + "_manta_out2.vcf"
			out3 = temp_dir + name + "_manta_out3.vcf"
			out4 = temp_dir + name + "_manta_out4.vcf"
			out5 = temp_dir + name + "_manta_out5.vcf"
			out6 = temp_dir + name + "_manta_out6.vcf"
			out7 = temp_dir + name + "_manta_out7.vcf"
			out8 = temp_dir + name + "_manta_out8.vcf"
			out9 = output_folder + name + "_manta_flagged.vcf"
			
			if parents != "FALSE":
				father_file = columns[3]
				mother_file = columns[4]
						
				father_manta = processed_dir + folder + "/structuralVariants/manta/" + father_file + "/results/variants/diploidSV.vcf.gz"
				mother_manta = processed_dir + folder + "/structuralVariants/manta/" + mother_file + "/results/variants/diploidSV.vcf.gz"
				
				
				print "#1 Filtering against father " + father_file
				run_flag_script(mode = "MANTA", input_file = input_manta, filter = father_manta, distance = dist, label = "PAT=1", output_file = out1)
				
				print "#2 Filtering against mother " + mother_file
				run_flag_script(mode = "MANTA", input_file = out1, filter = mother_manta, distance = dist, label = "MAT=1", output_file = out2)
				
			if parents == "FALSE":
				out2 = input_manta

			print "#3 Filtering against DELLY: " + delly_file
			run_flag_script(mode = "DELLY", input_file = out2, filter = delly_file, distance = dist, label = "DELLY=1", output_file = out3)
						
			print "#4 Filtering for ALT COV > 5"
			run_flag_script(filter_script = AC_script, mode = "AC", input_file = out3, threshold = "5", filter = "", distance = "", label = "ALTCOV", output_file = out4)
			
			print "#5 Filtering for ALT COV > 10"
			run_flag_script(filter_script = AC_script, mode = "AC", input_file = out4, threshold = "10", filter = "", distance = "", label = "ALTCOV", output_file = out5)

			print "#6 Filtering against 1000G"
			run_flag_script(mode = "1000G", input_file = out5, filter = VCF_1000G, distance = dist, label = "1000G=1", output_file = out6)
			
			print "#7 Filtering against DGV"
			run_flag_script(mode = "BED", input_file = out6, filter = VCF_DGV, distance = dist, label = "DGV=1", output_file = out7)

			print "#8 Filtering against VCF-explorer"
			run_flag_script(mode = "DB", input_file = out7, filter = manta_DB_filtered, distance = "0", label = "NOT_DB=1", output_file = out8)
						
			print "#9 Filtering against SpeedSeq Exclusion list"
			run_flag_script(mode = "BED", input_file = out8, filter = exclusion_file, distance = dist, label = "EXCL=1", output_file = out9)
			
f1.close()


