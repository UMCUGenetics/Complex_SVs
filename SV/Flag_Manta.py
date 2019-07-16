#!/usr/bin/python

import sys
import getopt
import re
import gzip

try:
    opts, args = getopt.getopt(sys.argv[1:], 'm:i:f:h:d:l:t:', ['mode=', 'input=', 'filter=', 'help=', 'distance=', 'label=', 'threshold='])
except getopt.GetoptError:
    usage()
    sys.exit(2)

for opt, arg in opts:
	if opt in ('-h', '--help'):
		usage()
		sys.exit(2)
	
	elif opt in ('-m', '--mode'):
		mode = arg
        
	elif opt in ('-i', '--input'):
		manta = arg
        
	elif opt in ('-f', '--filter'):
		filter_file = arg
    
	elif opt in ('-d', '--distance'):
		dist = int(arg)
		
	elif opt in ('-l', '--label'):
		label = arg
		
	elif opt in ('-t', '--threshold'):
		threshold = arg
		
	else:
		usage()
		sys.exit(2)

SVlist = dict()
mantaBND = dict()


if mode == "DELLY":
	with open(filter_file) as f1:
		for line in f1:
			line = line.rstrip()		
			if not line.startswith("#"):
				columns = line.split("\t")
				chr1 = columns[0]
				pos1 = int(columns[1])
				pos1a = pos1
				pos1b = pos1
			
				chr2 = re.match(".*CHR2=(\w+);", columns[7]).group(1)
				pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))			
				pos2a = pos2
				pos2b = pos2
			
				if "CIPOS" in columns[7]:
					cipos = map(int, re.match(".*CIPOS=-*(\d+,\d+)", columns[7]).group(1).split(","))
					pos1a -= cipos[0]
					pos1b += cipos[1]
				if "CIEND" in columns[7]:
					ciend = map(int, re.match(".*CIEND=-*(\d+,\d+)", columns[7]).group(1).split(","))
					pos2a -= ciend[0]
					pos2b += ciend[1]
			
				ori = "XX"
				if "CT=3to5" in columns[7]:
					ori = "TH"
				elif "CT=5to3" in columns[7]:
					ori = "HT"
				elif "CT=3to3" in columns[7]:
					ori = "TT"
				elif "CT=5to5" in columns[7]:
					ori = "HH"

				if chr2 < chr1:
					chr1, pos1a, pos1b, chr2, pos2a, pos2b = chr2, pos2a, pos2b, chr1, pos1a, pos1b
					ori = ori[::-1]
				elif chr2 == chr1 and pos2 < pos1:
					chr1, pos1a, pos1b, chr2, pos2a, pos2b = chr2, pos2a, pos2b, chr1, pos1a, pos1b
					ori = ori[::-1]
			
				
				if not "\t".join([chr1,chr2,ori]) in SVlist:
					SVlist["\t".join([chr1,chr2,ori])] = dict()
				if not pos1a in SVlist["\t".join([chr1,chr2,ori])]:
					SVlist["\t".join([chr1,chr2,ori])][pos1a] = dict()
				if not pos1b in SVlist["\t".join([chr1,chr2,ori])][pos1a]:
					SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b] = dict()
				if not pos2a in SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b]:
					SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b][pos2a] = dict()
				SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b][pos2a][pos2b] = 1

elif mode == "1000G":
	with open(filter_file) as f1:
		for line in f1:
			line = line.rstrip()		
			if not line.startswith("#"):
				columns = line.split(" ")
				chr1 = columns[0]
				pos1 = int(columns[1])
				pos1a = pos1
				pos1b = pos1
			
				#1000G SVs do not contain interchromosomal translocations, and does not contain a "CHR2=" field. Therefore chr2 should be the same as chr1
				#chr2 = re.match(".*CHR2=(\w+);", columns[8]).group(1)
				chr2 = chr1 
			
				# Not all SV calls contain an "END" field in the 1000G file.
				if "END" in columns[7]:
					pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))
				else:
					svlen = int(re.match(".*SVLEN=(\d+);", columns[7]).group(1))
					pos2 = pos1 + svlen
			
				pos2a = pos2
				pos2b = pos2
			
				if "CIPOS" in columns[7]:
					cipos = map(int, re.match(".*CIPOS=-*(\d+,\d+)", columns[7]).group(1).split(","))
					pos1a -= cipos[0]
					pos1b += cipos[1]
				if "CIEND" in columns[8]:
					ciend = map(int, re.match(".*CIEND=-*(\d+,\d+)", columns[7]).group(1).split(","))
					pos2a -= ciend[0]
					pos2b += ciend[1]
						
				if not "\t".join([chr1,chr2]) in SVlist:
					SVlist["\t".join([chr1,chr2])] = dict()
				if not pos1a in SVlist["\t".join([chr1,chr2])]:
					SVlist["\t".join([chr1,chr2])][pos1a] = dict()
				if not pos1b in SVlist["\t".join([chr1,chr2])][pos1a]:
					SVlist["\t".join([chr1,chr2])][pos1a][pos1b] = dict()
				if not pos2a in SVlist["\t".join([chr1,chr2])][pos1a][pos1b]:
					SVlist["\t".join([chr1,chr2])][pos1a][pos1b][pos2a] = dict()
				SVlist["\t".join([chr1,chr2])][pos1a][pos1b][pos2a][pos2b] = 1


elif mode == "MANTA":
	if re.search(".gz", filter_file):		
		f1 = gzip.open(filter_file)
	else:
		f1 = open(filter_file)
	
	for line in f1:
		line = line.rstrip()
		if not line.startswith("#"):
			columns = line.split("\t")
			chr1 = columns[0]
			pos1 = int(columns[1])
			chr2 = chr1
			pos2 = pos1
			pos1a = pos1
			pos1b = pos1
			ori = "XX"
			
			if "CIPOS" in columns[7]:
				cipos = map(int, re.match(".*CIPOS=-*(\d+,\d+)", columns[7]).group(1).split(","))
				pos1a -= cipos[0]
				pos1b += cipos[1]
			
			if "SVTYPE=DEL" in columns[7]:
				ori = "TH"
				pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))								
			elif "SVTYPE=INS" in columns[7]:
				ori = "TH"
				pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))				
			elif "SVTYPE=DUP" in columns[7]:
				ori = "HT"
				pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))
			elif "SVTYPE=INV" in columns[7]:
				pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))
				if "INV3" in columns[7]:
					ori = "TT"
				elif "INV5" in columns[7]:
					ori = "HH"
				else:
					ori = "YY"
			elif "SVTYPE=BND" in columns[7]:
				if not columns[2][:-1] in mantaBND:
					mantaBND[columns[2][:-1]] = [pos1a, pos1b]
				if columns[4].startswith("]"):
					ori = "HT"
					m = re.match("](\w+):(\d+)", columns[4])
					chr2 = m.group(1)
					pos2 = m.group(2)
				elif columns[4].startswith("["):
					ori = "HH"
					m = re.match("\[(\w+):(\d+)", columns[4])
					chr2 = m.group(1)
					pos2 = m.group(2)
				elif columns[4].endswith("["):
					ori = "TH"
					m = re.match(".*\[(\w+):(\d+)", columns[4])
					chr2 = m.group(1)
					pos2 = m.group(2)
				elif columns[4].endswith("]"):
					ori = "TT"
					m = re.match(".*](\w+):(\d+)", columns[4])
					chr2 = m.group(1)
					pos2 = m.group(2)
				else:
					ori = "ZZ"
			
			pos2a = pos2
			pos2b = pos2
			if columns[2][:-1] in mantaBND:
				pos2a = mantaBND[columns[2][:-1]][0]
				pos2b = mantaBND[columns[2][:-1]][1]
			if "CIEND" in columns[7]:
				ciend = map(int, re.match(".*CIEND=-*(\d+,\d+)", columns[7]).group(1).split(","))
				pos2a -= ciend[0]
				pos2b += ciend[1]
			if chr2 < chr1:
				chr1, pos1a, pos1b, chr2, pos2a, pos2b = chr2, pos2a, pos2b, chr1, pos1a, pos1b
				ori = ori[::-1]
			elif chr2 == chr1 and pos2 < pos1:
				chr1, pos1a, pos1b, chr2, pos2a, pos2b = chr2, pos2a, pos2b, chr1, pos1a, pos1b
				ori = ori[::-1]
			
			if not "\t".join([chr1,chr2,ori]) in SVlist:
				SVlist["\t".join([chr1,chr2,ori])] = dict()
			if not pos1a in SVlist["\t".join([chr1,chr2,ori])]:
				SVlist["\t".join([chr1,chr2,ori])][pos1a] = dict()
			if not pos1b in SVlist["\t".join([chr1,chr2,ori])][pos1a]:
				SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b] = dict()
			if not pos2a in SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b]:
				SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b][pos2a] = dict()
			SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b][pos2a][pos2b] = 1

elif mode == "DB":
	if re.search(".gz", filter_file):		
		f1 = gzip.open(filter_file)
	else:
		f1 = open(filter_file)
	
	for line in f1:
		line = line.rstrip()
		if not line.startswith("#"):
			columns = line.split("\t")
			if not re.search("DB1", columns[6]):
				chr1 = columns[0]
				pos1 = int(columns[1])
				chr2 = chr1
				pos2 = pos1
				pos1a = pos1
				pos1b = pos1
				ori = "XX"
		
				if "CIPOS" in columns[7]:
					cipos = map(int, re.match(".*CIPOS=-*(\d+,\d+)", columns[7]).group(1).split(","))
					pos1a -= cipos[0]
					pos1b += cipos[1]
		
				if "SVTYPE=DEL" in columns[7]:
					ori = "TH"
					pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))								
				elif "SVTYPE=INS" in columns[7]:
					ori = "TH"
					pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))				
				elif "SVTYPE=DUP" in columns[7]:
					ori = "HT"
					pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))
				elif "SVTYPE=INV" in columns[7]:
					pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))
					if "INV3" in columns[7]:
						ori = "TT"
					elif "INV5" in columns[7]:
						ori = "HH"
					else:
						ori = "YY"
				elif "SVTYPE=BND" in columns[7]:
					if not columns[2][:-1] in mantaBND:
						mantaBND[columns[2][:-1]] = [pos1a, pos1b]
					if columns[4].startswith("]"):
						ori = "HT"
						m = re.match("](\w+):(\d+)", columns[4])
						chr2 = m.group(1)
						pos2 = m.group(2)
					elif columns[4].startswith("["):
						ori = "HH"
						m = re.match("\[(\w+):(\d+)", columns[4])
						chr2 = m.group(1)
						pos2 = m.group(2)
					elif columns[4].endswith("["):
						ori = "TH"
						m = re.match(".*\[(\w+):(\d+)", columns[4])
						chr2 = m.group(1)
						pos2 = m.group(2)
					elif columns[4].endswith("]"):
						ori = "TT"
						m = re.match(".*](\w+):(\d+)", columns[4])
						chr2 = m.group(1)
						pos2 = m.group(2)
					else:
						ori = "ZZ"

				pos2a = pos2
				pos2b = pos2
				if columns[2][:-1] in mantaBND:
					pos2a = mantaBND[columns[2][:-1]][0]
					pos2b = mantaBND[columns[2][:-1]][1]
				if "CIEND" in columns[7]:
					ciend = map(int, re.match(".*CIEND=-*(\d+,\d+)", columns[7]).group(1).split(","))
					pos2a -= ciend[0]
					pos2b += ciend[1]
				if chr2 < chr1:
					chr1, pos1a, pos1b, chr2, pos2a, pos2b = chr2, pos2a, pos2b, chr1, pos1a, pos1b
					ori = ori[::-1]
				elif chr2 == chr1 and pos2 < pos1:
					chr1, pos1a, pos1b, chr2, pos2a, pos2b = chr2, pos2a, pos2b, chr1, pos1a, pos1b
					ori = ori[::-1]
		
				if not "\t".join([chr1,chr2,ori]) in SVlist:
					SVlist["\t".join([chr1,chr2,ori])] = dict()
				if not pos1a in SVlist["\t".join([chr1,chr2,ori])]:
					SVlist["\t".join([chr1,chr2,ori])][pos1a] = dict()
				if not pos1b in SVlist["\t".join([chr1,chr2,ori])][pos1a]:
					SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b] = dict()
				if not pos2a in SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b]:
					SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b][pos2a] = dict()
				SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b][pos2a][pos2b] = 1


elif mode == "TP":
	with open(filter_file) as f1:
		for line in f1:
			line = line.rstrip()		
			if not line.startswith("#"):
			
				columns = line.split("\t")
				chr1 = columns[1]
				pos1 = int(columns[2])
				pos1a = pos1
				pos1b = pos1
			
				chr2 = columns[3]
				pos2 = int(columns[4])
				pos2a = pos2
				pos2b = pos2

				if "+-" in columns[5]:
					ori = "TH"
				elif "-+" in columns[5]:
					ori = "HT"
				elif "++" in columns[5]:
					ori = "TT"
				elif "--" in columns[5]:
					ori = "HH"		
				else:
					ori = "XX"	
			
				if chr2 < chr1:
					chr1, pos1a, pos1b, chr2, pos2a, pos2b = chr2, pos2a, pos2b, chr1, pos1a, pos1b
					ori = ori[::-1]
				elif chr2 == chr1 and pos2 < pos1:
					chr1, pos1a, pos1b, chr2, pos2a, pos2b = chr2, pos2a, pos2b, chr1, pos1a, pos1b
					ori = ori[::-1]
				
				if not "\t".join([chr1,chr2,ori]) in SVlist:

					SVlist["\t".join([chr1,chr2,ori])] = dict()
				if not pos1a in SVlist["\t".join([chr1,chr2,ori])]:
					SVlist["\t".join([chr1,chr2,ori])][pos1a] = dict()
				if not pos1b in SVlist["\t".join([chr1,chr2,ori])][pos1a]:
					SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b] = dict()
				if not pos2a in SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b]:
					SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b][pos2a] = dict()
				SVlist["\t".join([chr1,chr2,ori])][pos1a][pos1b][pos2a][pos2b] = 1



elif mode == "BED":
	with open(filter_file) as f1:
		for line in f1:
			line = line.rstrip()
			if not line.startswith("#"):
				columns = line.split("\t")
				chr1 = columns[0]
				pos1 = int(columns[1])
				pos1a = pos1
				pos1b = pos1	

				chr2 = chr1 	
				pos2 = int(columns[2])
				pos2a = pos2
				pos2b = pos2
			
				if not "\t".join([chr1,chr2]) in SVlist:
					SVlist["\t".join([chr1,chr2])] = dict()
				if not pos1a in SVlist["\t".join([chr1,chr2])]:
					SVlist["\t".join([chr1,chr2])][pos1a] = dict()
				if not pos1b in SVlist["\t".join([chr1,chr2])][pos1a]:
					SVlist["\t".join([chr1,chr2])][pos1a][pos1b] = dict()
				if not pos2a in SVlist["\t".join([chr1,chr2])][pos1a][pos1b]:
					SVlist["\t".join([chr1,chr2])][pos1a][pos1b][pos2a] = dict()
				SVlist["\t".join([chr1,chr2])][pos1a][pos1b][pos2a][pos2b] = 1

else:
	print "#Mode not available. Available modes: DELLY, 1000G, MANTA, BED"


if re.search(".gz", manta):		
	f2 = gzip.open(manta)
else:
	f2 = open(manta)
	
for line in f2:
	line = line.rstrip()
	if not line.startswith("#"):
		columns = line.split("\t")
		chr1 = columns[0]
		pos1 = int(columns[1])
		ID = columns[2]
		REF = columns[3]
		ALT = columns[4]
		INFO = columns[7]
		chr2 = chr1
		pos2 = pos1
		pos1a = pos1
		pos1b = pos1
		ori = "XX"
		
		if "CIPOS" in columns[7]:
			cipos = map(int, re.match(".*CIPOS=-*(\d+,\d+)", columns[7]).group(1).split(","))
			pos1a -= cipos[0]
			pos1b += cipos[1]
		
		if "SVTYPE=DEL" in columns[7]:
			ori = "TH"
			pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))								
		elif "SVTYPE=INS" in columns[7]:
			ori = "TH"
			pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))				
		elif "SVTYPE=DUP" in columns[7]:
			ori = "HT"
			pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))
		elif "SVTYPE=INV" in columns[7]:
			pos2 = int(re.match(".*END=(\d+);", columns[7]).group(1))
			if "INV3" in columns[7]:
				ori = "TT"
			elif "INV5" in columns[7]:
				ori = "HH"
			else:
				ori = "YY"
		elif "SVTYPE=BND" in columns[7]:
			if not columns[2][:-1] in mantaBND:
				mantaBND[columns[2][:-1]] = [pos1a, pos1b]
			if columns[4].startswith("]"):
				ori = "HT"
				m = re.match("](\w+):(\d+)", columns[4])
				chr2 = m.group(1)
				pos2 = m.group(2)

			elif columns[4].startswith("["):
				ori = "HH"
				m = re.match("\[(\w+):(\d+)", columns[4])
				chr2 = m.group(1)
				pos2 = m.group(2)
			elif columns[4].endswith("["):
				ori = "TH"
				m = re.match(".*\[(\w+):(\d+)", columns[4])
				chr2 = m.group(1)
				pos2 = m.group(2)
			elif columns[4].endswith("]"):
				ori = "TT"
				m = re.match(".*](\w+):(\d+)", columns[4])
				chr2 = m.group(1)
				pos2 = m.group(2)
			else:
				ori = "ZZ"
		pos2a = pos2
		pos2b = pos2
		if columns[2][:-1] in mantaBND:
			pos2a = mantaBND[columns[2][:-1]][0]
			pos2b = mantaBND[columns[2][:-1]][1]
		if "CIEND" in columns[7]:
			ciend = map(int, re.match(".*CIEND=-*(\d+,\d+)", columns[7]).group(1).split(","))
			pos2a -= ciend[0]
			pos2b += ciend[1]
		if chr2 < chr1:
			chr1, pos1a, pos1b, chr2, pos2a, pos2b = chr2, pos2a, pos2b, chr1, pos1a, pos1b
			ori = ori[::-1]
		elif chr2 == chr1 and pos2 < pos1:
			chr1, pos1a, pos1b, chr2, pos2a, pos2b = chr2, pos2a, pos2b, chr1, pos1a, pos1b
			ori = ori[::-1]

		pos1a = pos1a - dist
		pos1b = pos1b + dist
		pos2a = pos2a - dist
		pos2b = pos2b + dist

		if mode == "1000G" or mode == "BED":
			if "\t".join([chr1,chr2]) in SVlist:
				for p1a in SVlist["\t".join([chr1,chr2])]:
					for p1b in SVlist["\t".join([chr1,chr2])][p1a]:
						for p2a in SVlist["\t".join([chr1,chr2])][p1a][p1b]:
							for p2b in SVlist["\t".join([chr1,chr2])][p1a][p1b][p2a]:
								if p1a <= pos1b and p1b >= pos1a and p2a <= pos2b and p2b >= pos2a:
									columns[7] = "".join([columns[7],";",label]) 
			print("\t".join([columns[i] for i in [0,1,2,3,4,5,6,7,8,9]]))


		else:
			if "\t".join([chr1,chr2,ori]) in SVlist:
				for p1a in SVlist["\t".join([chr1,chr2,ori])]:
					for p1b in SVlist["\t".join([chr1,chr2,ori])][p1a]:
						for p2a in SVlist["\t".join([chr1,chr2,ori])][p1a][p1b]:
							for p2b in SVlist["\t".join([chr1,chr2,ori])][p1a][p1b][p2a]:
								if p1a <= pos1b and p1b >= pos1a and p2a <= pos2b and p2b >= pos2a:
									columns[7] = "".join([columns[7],";",label]) 
								
			print("\t".join([columns[i] for i in [0,1,2,3,4,5,6,7,8,9]]))

	else:	
		print line

f1.close()
f2.close()



