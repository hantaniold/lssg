#benchmark.py
# 
# Sean Hogan
# Runs FPC on different datasets, with different table sizes.
# Measures runtime, maybe throughput as well, 
# for both decompression and compression.

import os 
import sys
import time
import math

def benchmark(filename,compfile,logfile,dataFolder):
# ct = compression time
# cl = compression level (table size in powers of 2 KB)
	cls = map(str,[1,3,5,7,9,11,13,15,17,19,21,23,25])
	for cl in cls:

# Dataset name
		ds = filename.lstrip(dataFolder+"/")

		print ds+" compressing, size "+cl+".....",
		cmd = "cat "+filename+" | bin/fpc "+cl+" > "+compfile
		it = time.time()
		os.system(cmd)
		ct = time.time() - it
		print ct
	
		print ds+" decompressing...",		
		cmd  = "cat "+compfile+" | bin/fpc > "+filename+".dc"
		it = time.time();	
		os.system(cmd)
		dt = time.time() - it
		print dt

		comp_size = float(os.path.getsize(filename)/math.pow(2,20))
		decomp_size = 	float(os.path.getsize(compfile)/math.pow(2,20))
		comp_tp	= comp_size / ct 
		decomp_tp = decomp_size / dt
		cr = comp_size/decomp_size
		print "comp ratio"+str(cr)
		logfile.write(\
			ds+"\t"\
			+"{:.4f}".format(ct)+\
			"\t"+"{:.2f}".format(comp_tp)+\
			"\t"+"{:.3f}".format(comp_size)+\
			"\t"+"{:.4f}".format(dt)+\
			"\t"+"{:.2f}".format(decomp_tp)+\
			"\t"+"{:.3f}".format(decomp_size)+\
			"\t"+cl+\
			"\t"+"{:.3f}".format(comp_size/decomp_size)+"\n"
		)
	
	

# Main
if (len(sys.argv) != 2) :
  print "Usage: python benchmark.py  <DATA>"
  exit(1)
dataFolder = sys.argv[1]

os.system("rm data/*.dc compressed-data/*")


file_list = os.listdir(dataFolder)
if not os.path.exists(dataFolder):
  print "Data folder "+dataFolder+" does not exist"
  exit(1)

if not os.path.exists("bin/fpc"):
  print "Expect fpc binary in directory bin."
  exit(1)

if not os.path.exists("fpc-results.log"):
	os.system("touch fpc-results.log")
logfile = open("fpc-results.log","a")
	

if not os.path.exists("compressed-data"):
	os.system("mkdir compressed-data")

logfile.write("------------------------\n")
logfile.write(time.ctime()+"\n")
logfile.write("------------------------\n")
logfile.write("NAME\t\t\t\tCT\t\tCTP\t\tDCS\t\tDT\t\tDTP\t\tCS\t\tTS\t\tCR\n")

for f in file_list:	
	benchmark("data/"+f,"compressed-data/"+f+".fpc",logfile,dataFolder)



