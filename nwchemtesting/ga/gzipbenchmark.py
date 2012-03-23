

import sys
import os
import time
import math


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
        
        
            



raw_input("Usage: inputfilelist inputdir logfile nrtrials\ny for ok\n")
inputfilelist = open(sys.argv[1],"r")
inputdir = sys.argv[2].rstrip("/")
logfile = open(sys.argv[3],"a")
nrTrials = int(sys.argv[4])



i = 0
j = 0
logfile.write("=============\n"+time.ctime()+"\n"+inputdir+"\n"+"==========="+"\ntype\ttr\tnr\tt\t\ttp\t\tuser\ttotal\n")
for f in inputfilelist:
    f = f.rstrip("\n") + ".bin"
    i += 1
    dcsize = os.path.getsize(inputdir+"/"+f)
    run_ctp = 0.0
    ctp = 0.0
    cr = 0.0
    run_ct = 0.0
    run_tt = 0.0
    total_ct = 0.0
    total_tt = 0.0
    while (j < nrTrials):
        j += 1
        os.system("( time gzip "+inputdir+"/"+f+" ) 2> gziptests.tmp")
        csize = float(os.path.getsize(inputdir+"/"+f+".gz"))
        os.system("gzip -d "+inputdir+"/"+f)  
        run_ct = parse_time("gziptests.tmp","compression")
        run_tt = parse_time("gziptests.tmp","total")
        run_ctp = ((dcsize/(run_ct)) / math.pow(2,20))

        cr = dcsize/csize
        if (j != 1):
            total_tt += run_tt
            total_ct += run_ct
            ctp += run_ctp
        else:
            print "Skippin first trial"
            
        logfile.write("trial\t "+str(j)+"\t"+str(i)+"\t"+"{:.4f}".format(cr)+"\t"+"{:.4f}".format(run_ctp)+"\t"+"{:.4f}".format(run_ct)+"\t""{:.4f}".format(run_tt)+"\n")
        print("Finished trial "+str(j)+" "+f+"\t cr: "+str(cr))
    j = 0
    
    logfile.write("av\t\t-\t"+str(i)+"\t"+"{:.4f}".format(cr)+"\t"+"{:.4f}".format(ctp/(nrTrials - 1))+"\t"+"{:.4f}".format(total_ct/(nrTrials - 1))+"\t"+"{:.4f}".format(total_tt/(nrTrials-1))+"\n")
    
    
