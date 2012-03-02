import os
import sys

#Takes a list of files changes them into binary representations

res = raw_input("The input should be a list of TCE input \n \
with 5 lines of info we don't need, then lines of [number] [some double].\n \
 The output will be all of the input filenames plus \".bin\" appended, \n \
 written as doubles for input into the GA program.\n \
 If this is okay type \"y\"\n")

if res != "y":
    exit(0)

filenames = open(sys.argv[1],"r")

for filename in filenames:
    # Remove the leading 5 lines
    filename = filename.rstrip('\n')
    os.system("tail -n +6 "+filename+" > "+filename+"mytmp")
    # Rename the temp file to the original file.
    # Convert into binary using ga2bin
    os.system("./ga2bin "+filename+"mytmp "+filename+".bin")
    os.system("rm "+filename+"mytmp")
