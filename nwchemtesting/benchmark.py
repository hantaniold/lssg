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

# Global
nrTrials = 1

def parse_time(f, tp):
    fi = open(f,"r")
    usertime = 0.0
    totaltime = 0.0
    user = 0.0
    sys = 0.0
    for line in fi:
        line = line.split()
        user = float(line[0].split("u")[0])
        sys = float(line[1].split("s")[0])
        if (tp == "compression"):
            return user
        else:
            return (user + sys)
 

def benchmark(filename,compfile,logfile,dataFolder,i):
# ct = compression time
# cl = compression level (table size in powers of 2 KB)
	cls = map(str,[15])
#   cls = map(str,[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25])
#	cls = map(str,[1,4,7,10,12,13,14,15,16,17,19,22,25])

	for cl in cls:

# Dataset name
		ds = filename.lstrip(dataFolder+"/")

        j = 0
        ct = 0  # user time
        tct = 0
        all_ct = 0  # user + sys
        all_tct = 0
        run_ctp = 0.0
        total_ctp = 0.0

        dt = 0
        decomp_tp = 0
        
        
        while (j < nrTrials):
            j += 1
            cmd = "(time cat "+filename+" | bin/fpc "+cl+" > "+compfile+") 2> fpcbench.tmp"
            os.system(cmd)

            ct = parse_time("fpcbench.tmp","compression")
            all_ct = parse_time("fpcbench.tmp","total")
            
        
            print ds+" decompressing...",		
            cmd  = "cat "+compfile+" | bin/fpc > "+filename+".dc"
            it = time.time();	
            os.system(cmd)
            dt += time.time() - it
            print dt
        
            comp_size = float(os.path.getsize(filename)/math.pow(2,20))
            run_ctp	= comp_size / ct 
            if (j != 1):
                total_ctp += run_ctp
                tct += ct
                all_tct += all_ct
            else:
                print "Skipping first run"
            decomp_size = 	float(os.path.getsize(compfile)/math.pow(2,20))
            decomp_tp += decomp_size / dt

            cr = comp_size/decomp_size
            logfile.write("trial\t "+str(j)+"\t"+str(i)+"\t"+"{:.4f}".format(cr)+"\t"+"{:.4f}".format(run_ctp)+"\t"+"{:.4f}".format(ct)+"\t""{:.4f}".format(tct)+"\n")
            print "trial "+str(j)+" comp ratio:\t"+str(cr)

        print("av tp: "+str(comp_tp/nrTrials))
        print("av ct: "+str(tct/nrTrials))
        logfile.write("av\t\t-\t"+str(i)+"\t"+"{:.4f}".format(cr)+"\t"+"{:.4f}".format(total_ctp/(nrTrials - 1))+"\t"+"{:.4f}".format(all_ct/(nrTrials - 1))+"\t"+"{:.4f}".format(all_tct/(nrTrials-1))+"\n")
	
	

# Main
if (len(sys.argv) != 5) :
  print "Usage: python benchmark.py  <DATA> <LOGFILE> <FILENAMELIST> <NRTRIALS>"
  exit(1)
dataFolder = sys.argv[1]

os.system("rm "+dataFolder+"/*.dc compressed-data/*")


if not os.path.exists(dataFolder):
  print "Data folder "+dataFolder+" does not exist"
  exit(1)
file_list = open(sys.argv[3],"r")

if not os.path.exists("bin/fpc"):
  print "Expect fpc binary in directory bin."
  exit(1)

if not os.path.exists("fpc-results.log"):
	os.system("touch fpc-results.log")
logfile = open(sys.argv[2],"a")
	

if not os.path.exists("compressed-data"):
	os.system("mkdir compressed-data")

logfile.write("=============\n"+time.ctime()+"\n"+dataFolder+"\n"+"==========="+"\ntype\ttr\tnr\tt\t\ttp\t\tuser\ttotal\n")

nrTrials = int(sys.argv[4])

i = 0
for f in file_list:	
    i += 1
    f = f.rstrip("\n") + ".bin"
    benchmark(dataFolder+"/"+f,"compressed-data/"+f+".fpc",logfile,dataFolder,i)

os.system("rm "+dataFolder+"/*.dc")
print("COMPLETE")


