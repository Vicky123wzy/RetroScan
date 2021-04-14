#!/usr/bin/env python
import sys
import argparse
import os
import subprocess
import re
from concurrent.futures import ThreadPoolExecutor
import time

USAGE = "\nusage: python %s tmp_path identity\n" % sys.argv[0]

if len(sys.argv) != 12:
    print(USAGE)
    sys.exit()

lossintronfile = sys.argv[1]
workdir = sys.argv[2]
kaks_single_bash = sys.argv[3]
genome_file = sys.argv[4]
protein_file = sys.argv[5]
bedtools = sys.argv[6]
kaksmethod = sys.argv[7]
seqkit = sys.argv[8]
diamond = sys.argv[9]
clustalw2 = sys.argv[10]
thread = sys.argv[11]


def callCMD(cmd):
    FNULL = open(os.devnull, 'w')
    subprocess.call(cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    FNULL.close()

def kaks_single_shell(kaks_single_bash,tmp_path,genome_file,protein_file,lossintronfile,cdsfile,bedtools,kaksmethod,seqkit,diamond,clustalw2):
    kaks_single_cmd = "bash " + kaks_single_bash + " " + tmp_path + " " +  genome_file + " " + protein_file + " " + \
                      lossintronfile + " " + cdsfile + " " + bedtools + " " + kaksmethod + " " + seqkit + " " + diamond + " " + clustalw2
    callCMD(kaks_single_cmd)
    time.sleep(0.001)


list_a =[]
with open(lossintronfile, 'r') as file:
        ls = file.readlines()[1:]
        for line in ls:
            line = line.rstrip("\n")
            array = line.split("\t")
            list_a.append(array[0])
            kaks_dir = "kaks_result/" + str(array[0])
            os.makedirs(os.path.join(workdir, kaks_dir))
            OUT = open(os.path.join(workdir, kaks_dir,"single_loss_intron.txt") , 'w')
            OUT.write(line)
            OUT.close()

with ThreadPoolExecutor(int(thread)) as executor:
    for i in list_a:
        kaks_dir = "kaks_result/" + str(i)
        working_path = os.path.join(workdir, kaks_dir)
        inputfile = os.path.join(working_path,"single_loss_intron.txt")
        cdsfile = os.path.join(workdir, "cds.fasta")
        executor.submit(kaks_single_shell,kaks_single_bash,working_path,genome_file,protein_file,inputfile,cdsfile,bedtools,kaksmethod,seqkit,diamond,clustalw2)



