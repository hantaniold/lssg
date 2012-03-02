#Takes a list of files and computes successive deltas
#Writes this output to some file (compression info) decompression info with FPC and GZIP.


import os
import sys



### MAIN ###

if (len(sys.argv) < 3):
    print("usage: [input file list] [input dir] [output dir] [output prefix] <-np #nrprocs> <-f cutoff exponent>");
    exit(0)

inputFileHandle = open(sys.argv[1],"r")
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
print "=================================="
ans = raw_input("enter y if this is okay\n>>> ")
if (ans != "y"):
  print("Exiting.")
  exit(0)


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
    os.system(mpiCmd+subt+minu+output)
    
# Compress the files.
# Decompress.
# Recover.

for i in range(len(inputFilesList) - 1):
    subt = "-s "+inputDir+inputFilesList[i]+".bin "
    delt = "-d "+outputDir+outputPrefix+"."+str(i+2)+"-"+str(i+1)+".del "
    output = "-o "+outputDir+inputFilesList[i+1]+".rec"
    os.system(mpiCmd+" -r "+subt+delt+output)
    print(mpiCmd+" -r "+subt+delt+output)

    
