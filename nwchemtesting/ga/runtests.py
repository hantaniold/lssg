#Takes a list of files and computes successive deltas
#Writes this output to some file (compression info) decompression info with FPC and GZIP.


import os
import sys
import time
import math



### MAIN ###

if (len(sys.argv) < 3):
    print("usage: [input file list] [input dir] [output dir] [output prefix] <-np #nrprocs> <-f cutoff exponent>");
    exit(0)

inputFileHandle = open(sys.argv[1],"r")
logfile = open(sys.argv[1]+".log","a")
inputFilesList = []
for line in inputFileHandle:
    inputFilesList += [line.rstrip("\n")]

inputDir = sys.argv[2].rstrip("/")+"/"
outputDir = sys.argv[3].rstrip("/")+"/"
outputPrefix = sys.argv[4]


np = 1
cutoff = ""
# Parse rest of the optional params
for i in range(4,len(sys.argv)):
    if (sys.argv[i] == "-np"):
        np = int(sys.argv[i+1])
    elif (sys.argv[i] == "-f"):
        cutoff = " -f "+sys.argv[i+1]+" "


print "=================================="
print "list of input files: "+sys.argv[1]
print "input dir: "+sys.argv[2]
print "output directory: "+sys.argv[3]
print "output prefix: "+sys.argv[4]
print "cutoff args: "+cutoff
print "nr procs: "+str(np)
print "logfile: "+sys.argv[1]+".log"
print "=================================="
ans = raw_input("enter y if this is okay\n>>> ")
if (ans != "y"):
  print("Exiting.")
  exit(0)

logfile.write("------------------------\n")
logfile.write(time.ctime()+"\n")
logfile.write("------------------------\n")
logfile.write("CUTOFF: "+cutoff+" NRPROCS: "+str(np)+"\n")
logfile.write("timestep\tCR\ttime\tthroughput\n")


#Compute deltas of 2-1, 3-2
#Compress/decompress these deltas and time and whatnot
#Compute recovery and do the times or something and check sizes


mpiCmd = "mpirun -np "+str(np)+" -machinefile machines.txt ./ga-delta.x "+cutoff

# Run delta compression. Write out files with successive
# numeric prefixes: outputPrefix.n-n-1.delta
for i in range(len(inputFilesList) - 1):
    
    subt = "-s "+inputDir+inputFilesList[i]+".bin "
    minu = "-m "+inputDir+inputFilesList[i+1]+".bin "
    output =  "-o "+outputDir+outputPrefix+"."+str(i+2)+"-"+str(i+1)+".del"
    print(mpiCmd+subt+minu+output)
    tic = time.time()
    os.system(mpiCmd+subt+minu+output)
    toc = time.time()
    comp_size = float(os.path.getsize(minu.split()[1]))
    decomp_size = float(os.path.getsize(output.split()[1]))
    comp_ratio = comp_size/decomp_size
    ct = toc - tic
    comp_tp = (comp_size/ct)/math.pow(2,20)
    logfile.write(str(i+1)+"\t"+"{:.4f}".format(comp_ratio)+\
        "\t"+"{:.4f}".format(ct)+"\t"+"{:.4f}".format(comp_tp)+"\n")
    

print("DONE")
logfile.close()
inputFileHandle.close()
exit(1)
    
# Recover.

for i in range(len(inputFilesList) - 1):
    subt = "-s "+inputDir+inputFilesList[i]+".bin "
    delt = "-d "+outputDir+outputPrefix+"."+str(i+2)+"-"+str(i+1)+".del "
    output = "-o "+outputDir+inputFilesList[i+1]+".rec"
    os.system(mpiCmd+" -r "+subt+delt+output)
    print(mpiCmd+" -r "+subt+delt+output)

    
