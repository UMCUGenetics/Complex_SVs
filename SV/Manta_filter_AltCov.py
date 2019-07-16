import sys
import re
import gzip
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:h:l:t:', ['input=', 'help=', 'label=', 'threshold='])
except getopt.GetoptError:
    usage()
    sys.exit(2)

for opt, arg in opts:
	if opt in ('-h', '--help'):
		usage()
		sys.exit(2)
	       
	elif opt in ('-i', '--input'):
		manta = arg
        
	elif opt in ('-l', '--label'):
		label = arg
		
	elif opt in ('-t', '--threshold'):
		threshold = arg
		
	else:
		usage()
		sys.exit(2)

with open(manta) as f1:
	for line in f1:
		line = line.rstrip()
		if not line.startswith("#"):
			columns = line.split("\t")
			
			genotype = columns[9]

			#print genotype
			
			PR = genotype.split(":")[4]
			ALTCOV = int(PR.split(",")[1])

			#print PR
			#print ALTCOV

			if ALTCOV > int(threshold):
				columns[7] = "".join([columns[7],";",label, ">",str(threshold)]) 
				print("\t".join([columns[i] for i in [0,1,2,3,4,5,6,7,8,9]]))
			else:
				print line
									


		#else:	
			#print line

f1.close()
