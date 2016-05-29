import commands
import shutil
import os
import sys
import time
import linecache

import math

do_submit = True

delete_help_file = True # if True, delete all current jobs_*, out, err files

Ncell = 2  # e.g. 5 supercell = 10x10 lattice
mus = [0.36]

Us = [-6]
Ds = [0.3]
tp = 0.4

seeds = [3234567]#, 3234567]#, 3234567, 4234567, 5234567, 6234567, 7234567, 8234567, 9234567]

dtaus   = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
Ls      = [20]
#Ls      = [20, 30, 40, 50, 60, 70, 80]#, 100, 120, 140]
ntry    = [1,4,4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
ntry2   = [1,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2]

tdm = 1
HSF = -1      # -1 = random, 1 = from file 
nbin = 10      # only applies for non-MPI run
FTphy = 0    # if do FT for phy0
FTtdm = 0
SelfE  = 0    # if self-energy
Dsqy   = 0    # from curr-curr correlation

hours   = [0, 24, 24, 24, 24, 24, 24, 10, 3, 3, 3, 5, 6, 5, 6, 7, 8]
minutes = [12, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0,0, 0, 0, 0]

CPUS = [1]  

def clean():

    for i in range(0,len(temp)):

        if os.path.exists("./T_" + str(temp[i])):
            shutil.rmtree("./T_" + str(temp[i]))

        cmd = "mkdir T_" + str(temp[i])
        os.system(cmd)

def prepare_input_file(filename, outputname, i, s, m):

    file = open(filename, "r")
    text = file.read()
    file.close()

    text = text.replace("OUTPUT"    , str(outputname))
    text = text.replace("GEOM"      , str(geomfile))

    if(Ls[i] >= 50):
        text = text.replace("NMEASval"    , str(2500))
    else:
        text = text.replace("NMEASval"    , str(2500))

    if(Ls[i] >= 50):
        text = text.replace("NWARMval"   , str(1000))
    else:
        text = text.replace("NWARMval"   , str(1000))

    text = text.replace("SEEDval"    , str(seeds[s]))
    text = text.replace("MUval"      , str(mus[m]))

    text = text.replace("Lval"    , str(Ls[i]))
    text = text.replace("DTAUval" , str(dtaus[i]))

    text = text.replace("NTRYval"   , str(ntry[i]))
    text = text.replace("NTRY2val"  , str(ntry2[i]))
    text = text.replace("TDMval"    , str(tdm))

    text = text.replace("NORTHval"  , str(10))

    text = text.replace("HSFval"    , str(HSF))
    text = text.replace("NBINval"   , str(nbin))

    text = text.replace("FTPHY" , str(FTphy))
    text = text.replace("SELFEval" , str(SelfE))
    text = text.replace("DSQYval" ,  str(Dsqy))

    file = open(filename, "w")
    file.write(text)
    file.close()

def prepare_geom_file(geomfile, iU, iD):

    file = open(geomfile, "r")
    text = file.read()
    file.close()

    text = text.replace("NCELL",  str(Ncell))
    text = text.replace("Uval" ,  str(Us[iU]))
    text = text.replace("Dval" ,  str(Ds[iD]))
    text = text.replace("tp"   ,  str(tp))

    file = open(geomfile, "w")
    file.write(text)
    file.close()

cmd = "rm jobs_*"
if(delete_help_file):
    os.system(cmd)

#cmd = "rm -r run_*"
#if(delete_help_file):
#    os.system(cmd)

#cmd = "rm out_*"
#if(False and delete_help_file):
#    os.system(cmd)

cmd = "rm error_*"
if(False and delete_help_file):
    os.system(cmd)

for iD in range(0, len(Ds)):
      for iU in range(0, len(Us)):
   
        for m in range(0, len(mus)):
 
	    lattice  = "U"+str(Us[iU])+"_D"+str(Ds[iD])+"_tp"+str(tp)+"_mu"+str(mus[m])+"_N"+str(8*Ncell**3)
       	    geomfile = "geomU"+str(Us[iU])+"_D"+str(Ds[iD])+"_tp"+str(tp)+"_N"+str(8*Ncell**3)
            print lattice  

            cmd = "cp g_template" + " "+ geomfile
            print cmd
            os.system(cmd)

            prepare_geom_file(geomfile, iU, iD)

            for i in range(0, len(Ls)):

                beta = Ls[i]*dtaus[i]
                print "L=", Ls[i], "dtau=", dtaus[i], "beta = ", beta

                dir = "./"+lattice+"_be" + str(beta)

                if not os.path.exists(dir):
                  cmd = "mkdir " + dir
                  os.system(cmd)
                  cmd = "cp ggeom " + dir
                  os.system(cmd)
                  cmd = "cp " + geomfile + " "+ dir
                  os.system(cmd)

                for s in range(0,len(seeds)):

                    input_file_name    = dir + "/input_seed" + str(seeds[s])
                    data_file_name     = lattice + "_be" + str(beta)+ "_s" + str(seeds[s])+ "_mpi" + str(CPUS[0]*8)

                    cmd = "cp ./in_template " + input_file_name + ";"
                    os.system(cmd)
                    prepare_input_file(input_file_name, data_file_name, i, s, m)

                    batch_str = ""
                    batch_str = batch_str + "aprun -B ./ggeom " + dir + "/input_seed"+ str(seeds[s]) +"\n"

                    file = open("batch_script.slm", "r")
                    text = file.read()
                    file.close()

                    jobs_file_name = dir + "/jobs_"+lattice+"_be"+ str(beta)+ "_s"+str(seeds[s])+".slm"
                
                    text = text.replace("CPUS"  , str(CPUS[0]))
                    text = text.replace("MPI"   , str(CPUS[0]*8))
                    text = text.replace("HOUR"  , str(hours[i]))
                    text = text.replace("MINUTE", str(minutes[i]))
                    text = text.replace("BETA" , str(beta))
                    text = text.replace("SEED" , str(seeds[s]))
                    text = text.replace("JOBS" , str(batch_str))
                    text = text.replace("LATTICE" , lattice)

                    if(not batch_str == ""):
                    	file = open(jobs_file_name, "w")
                    	file.write(text)
                    	file.close()
     
                    if(do_submit):
                    	cmd = "sbatch " + jobs_file_name
                     	os.system(cmd)

