#!/usr/bin/env python
import sys

USAGE = "\nusage: python %s blast retropro retrodna outdna outpro description retroprotein\n" % sys.argv[0]

if len(sys.argv) != 8:
    print(USAGE)
    sys.exit()

score = 0
with open(sys.argv[1], 'r') as f:
    for line in f:
        line = line.rstrip("\n")
        array = line.split("\t")
        if float(array[11]) > float(score) :
            score = array[11]
            retro_name = array[0]

retro_proseq = ""
with open(sys.argv[2], 'r') as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            retro_pro = line.replace(">","").replace(" ", "")
        elif retro_pro == retro_name:
                retro_proseq += line

output_file = open(sys.argv[7], "w")
output_file.write(">" + retro_pro + "\n" + retro_proseq + "\n")
output_file.close()


retro_dnaseq = ""
with open(sys.argv[3], 'r') as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            retro = line.replace(">","").replace(" ", "")
        else:
            retro_dnaseq += line

def getStrInfo(str,target):
    str_site = []
    for index,value in enumerate(str):
        if target == value:
            str_site.append(index)
    return str_site

frame_site = retro_name.split("=")[1]
list_dna = list(retro_dnaseq.strip())
complement = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A', 'a': 'T', 'g': 'C', 'c': 'G', 't': 'A', 'N': 'N', 'n': 'N'}
rev_seq = ''
strand = ""
des = ""
if frame_site == "1":
    strand = "+"
elif frame_site == "2":
    strand = "+"
    del list_dna[0]
elif frame_site == "3":
    strand = "+"
    del list_dna[0:2]
elif frame_site == "-1":
    strand = "-"
    for i in list_dna:
        rev_seq += complement[i]
    rev_seq = rev_seq[::-1]
    list_dna = list(rev_seq)
elif frame_site == "-2":
    strand = "-"
    for i in list_dna:
        rev_seq += complement[i]
    rev_seq = rev_seq[::-1]
    list_dna = list(rev_seq)
    del list_dna[0]
elif frame_site == "-3":
    strand = "-"
    for i in list_dna:
        rev_seq += complement[i]
    rev_seq = rev_seq[::-1]
    list_dna = list(rev_seq)
    del list_dna[0:2]

if '*' not in retro_proseq:
    des = "intact"
    outdna = "".join(str(i) for i in list_dna)
    outpro = retro_proseq

elif '*' in retro_proseq:
    des = "pseudogene"
    stop_site = getStrInfo(retro_proseq,"*")
    outpro = retro_proseq.replace("*","")
    remove_site = []
    for i, val in enumerate(stop_site):
        remove_site.append(int(val) * 3)
        remove_site.append(int(val) * 3 + 1)
        remove_site.append(int(val) * 3 + 2)

    list_dna = [list_dna[i] for i in range(0, len(list_dna), 1) if i not in remove_site]
    outdna = "".join(str(i) for i in list_dna)

output_file = open(sys.argv[4], "w")
output_file.write(">" + retro + "\n" + outdna + "\n")
output_file.close()

output_file = open(sys.argv[5], "w")
output_file.write(">" + retro + "\n" + outpro + "\n")
output_file.close()

output_file = open(sys.argv[6], "w")
output_file.write(retro + "\t" + strand + "\t" + des + "\n")
output_file.close()
