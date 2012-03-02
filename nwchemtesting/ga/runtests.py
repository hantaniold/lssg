#Takes a list of files and computes successive deltas
#Writes this output to some file (compression info) decompression info with FPC and GZIP.


import os
import sys



### MAIN ###

inputFileHandle = open(sys.argv[1],"r")
inputFilesList = []
for line in inputFileHandle:
    inputFilesList += [line.rstrip("\n")]

#Compute deltas of 2-1, 3-2
#Compress/decompress these deltas and time and whatnot
#Compute recovery and do the times or something and check sizes
