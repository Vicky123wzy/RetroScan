#!/usr/bin/env python
from __future__ import print_function
import sys
import argparse
import os
import subprocess
import re
import zipfile
import tarfile
import logging
from logging import handlers
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

# log setting
class Logger(object):
    level_relations = {
        'debug': logging.DEBUG,
        'info': logging.INFO,
        'warning': logging.WARNING,
        'error': logging.ERROR,
        'crit': logging.CRITICAL
    }

    def __init__(self, filename, level='debug', when='D', backCount=3,
                 fmt='[%(asctime)s] - %(levelname)s: %(message)s'):
        self.logger = logging.getLogger(filename)
        format_str = logging.Formatter(fmt, datefmt='%Y-%m-%d %H:%M:%S')
        self.logger.setLevel(self.level_relations.get(level))
        sh = logging.StreamHandler()
        sh.setFormatter(format_str)
        th = handlers.TimedRotatingFileHandler(filename=filename, when=when, backupCount=backCount, encoding='utf-8')
        th.setFormatter(format_str)
        self.logger.addHandler(sh)
        self.logger.addHandler(th)

# log = Logger('all.log',level='debug')
# log.logger.debug('debug')
# log.logger.info('info')
# log.logger.warning('warning')
# log.logger.error('error')
# log.logger.critical('critical')


def callCMD(cmd):
    FNULL = open(os.devnull, 'w')
    subprocess.call(cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    FNULL.close()


# check input files and working dir
def check_input_data(genome_file, gff_file, working_dir):
    # creat working_dir
    working_dir = os.path.abspath(working_dir)
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
        if not os.path.exists(working_dir):
            print("ERROR: It is not possible to create result directory.")
            sys.exit(0)

    if not working_dir.endswith('/'):
        working_dir += '/'

    if not os.path.exists(working_dir + '/retroscan_results'):
        os.makedirs(working_dir + '/retroscan_results')
 

    # check if genome file and gff file exist
    genome_file = os.path.abspath(genome_file)
    if not os.path.exists(genome_file):
        print('ERROR: The genome file does not exist.')
        sys.exit(0)
    elif genome_file.endswith(".fasta") or genome_file.endswith(".fa") \
            or genome_file.endswith(".fas") or genome_file.endswith(".faa") or genome_file.endswith(".fna"):
        genome_file = os.path.abspath(genome_file)
    elif genome_file.endswith(".tar.gz"):
        t = tarfile.open(genome_file)
        tmp_path = os.path.join(working_dir, "results/tmp")
        t.extractall(path=tmp_path)
        genome_file = os.path.join(tmp_path, t.getnames()[0])
        t.close()
    elif genome_file.endswith(".zip"):
        zip_file = zipfile.ZipFile(genome_file)
        zip_list = zip_file.namelist()[0]
        tmp_path = os.path.join(working_dir, "results/tmp")
        zip_file.extract(zip_list, tmp_path)
        genome_file = os.path.join(tmp_path, zip_list)
        zip_file.close()
    else:
        print("ERROR: The " + genome_file + " is not right!")
        sys.exit(0)

    gff_file = os.path.abspath(gff_file)
    if not os.path.exists(gff_file):
        print('ERROR: The gff file does not exist.')
        sys.exit(0)
    elif gff_file.endswith(".gff") or gff_file.endswith(".gff3"):
        gff_file = os.path.abspath(gff_file)
    elif gff_file.endswith(".tar.gz"):
        t = tarfile.open(gff_file)
        tmp_path = os.path.join(working_dir, "results/tmp")
        t.extractall(path=tmp_path)
        gff_file = os.path.join(tmp_path, t.getnames()[0])
        t.close()
    elif gff_file.endswith(".zip"):
        zip_file = zipfile.ZipFile(gff_file)
        zip_list = zip_file.namelist()[0]
        tmp_path = os.path.join(working_dir, "results/tmp")
        zip_file.extract(zip_list, tmp_path)
        gff_file = os.path.join(tmp_path, zip_list)
        zip_file.close()
    else:
        print("ERROR: The " + gff_file + " is not right!")
        sys.exit(0)

    return genome_file, gff_file, working_dir


# calculate length of proteins
def PRO_SEQ(pro_file):
    pro_seq = {}
    pro_len = {}
    for line in pro_file:
        line = line.rstrip("\n")
        if line.startswith('>'):
            mRNA = line.split()[0].replace(">", "")
            seq = ''
        else:
            seq += line
            pro_seq[mRNA] = seq
    for key in pro_seq:
        pro_len[key] = len(pro_seq[key])
    return pro_seq, pro_len


# extract longest protein sequence for each gene which has three or more exons
def gene_trans_info(tmp_path, protein_path, gff_file):
    pro_file = open(protein_path, "r")
    pro_seq, pro_len = PRO_SEQ(pro_file)

    # file of gene_trans
    mRNA_out = {}
    cds_site = {}
    exon_site = {}
    gene_trans = open(os.path.join(tmp_path, "gene_trans_info.txt"), "w")
    gene_trans.write(
        "gene\tmRNA\tmRNA_chr\tmRNA_strand\tmRNA_start\tmRNA_end\tprotein_length\tcds_num\tcds_site\texon_site\tpep_site\n")
    with open(gff_file, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            array = line.split("\t")
            if (len(array) > 7 and array[2] == "mRNA") or (len(array) > 7 and array[2] == "transcript"):
                des = re.split('[;,]', array[8])
                for j in range(len(des)):
                    if "Parent=" in des[j] or "parent=" in des[j]:
                        gene = des[j].split('=')[-1]
                mRNA = des[0].split('=')[-1]
                if mRNA in pro_len and mRNA not in mRNA_out:
                    mRNA_out.setdefault(mRNA, []).append(gene)  # gene
                    mRNA_out.setdefault(mRNA, []).append(array[0])  # chr
                    mRNA_out.setdefault(mRNA, []).append(array[6])  # mRNA_strand
                    mRNA_out.setdefault(mRNA, []).append(array[3])  # mRNA_start
                    mRNA_out.setdefault(mRNA, []).append(array[4])  # mRNA_end
            elif len(array) > 7 and array[2] == "CDS":
                des = re.split('[;,]', array[8])
                for j in range(len(des)):
                    if "Parent=" in des[j] or "parent=" in des[j]:
                        cds_parent = des[j].split('=')[-1]
                if cds_parent not in cds_site:
                    cds_site.setdefault(cds_parent, []).append(int(1))  # cds_num
                    cds_site.setdefault(cds_parent, []).append(array[3] + "-" + array[4])  # cds_site
                elif cds_parent in cds_site:
                    cds_site[cds_parent][0] = cds_site[cds_parent][0] + 1  # cds_num
                    cds_site[cds_parent][1] = cds_site[cds_parent][1] + "," + array[3] + "-" + array[4]  # cds_site
            elif len(array) > 7 and array[2] == "exon":
                des = re.split('[;,]', array[8])
                for j in range(len(des)):
                    if "Parent=" in des[j] or "parent=" in des[j]:
                        exon_parent = des[j].split('=')[-1]
                if exon_parent not in exon_site:
                    exon_site[exon_parent] = array[3] + "-" + array[4]
                elif exon_parent in exon_site:
                    exon_site[exon_parent] = exon_site[exon_parent] + "," + array[3] + "-" + array[4]

    if not exon_site:
        for key in cds_site:
            exon_site[key] = cds_site[key][1]

    for key in pro_seq:
        if key in mRNA_out and key in cds_site and key in exon_site:
            pep_site = []
            pep_start = 0
            remainder = 0
            cds_site_list = cds_site[key][1].split(",")
            if len(cds_site_list) == 1:
                pep_site.append("NA")
            elif len(cds_site_list) > 1:
                if mRNA_out[key][2] == "-":
                    cds_site_list.sort(reverse=True)
                elif mRNA_out[key][2] == "+":
                    cds_site_list.sort()
                for num in cds_site_list:
                    len_cds = int(num.split("-")[1]) - int(num.split("-")[0]) + 1
                    pep_start = (len_cds + int(remainder)) // 3 + pep_start
                    remainder = (len_cds + int(remainder)) % 3
                    pep_site.append(pep_start)
                pep_site.pop()
            gene_trans.write(str(mRNA_out[key][0]) + "\t" + str(key) + "\t" + str(mRNA_out[key][1]) + "\t" +
                             str(mRNA_out[key][2]) + "\t" + str(mRNA_out[key][3]) + "\t" + str(mRNA_out[key][4]) + "\t" +
                             str(pro_len[key]) + "\t" + str(cds_site[key][0]) + "\t" + str(cds_site[key][1]) + "\t" +
                             str(exon_site[key]) + "\t" + ",".join(str(i) for i in pep_site) + "\n")
    gene_trans.close()

    # extract longest protein which has three or more exons
    longest_len = {}
    longest_mRNA = {}
    with open(os.path.join(tmp_path, "gene_trans_info.txt"), 'r') as f:
        lines = f.readlines()[1:]
        for line in lines:
            line = line.rstrip("\n")
            array = line.split("\t")
            if int(array[7]) > 2:
                if array[0] not in longest_len:
                    longest_len[array[0]] = array[6]
                    longest_mRNA[array[0]] = array[1]
                elif array[0] in longest_len and int(array[6]) > int(longest_len[array[0]]):
                    longest_len[array[0]] = array[6]
                    longest_mRNA[array[0]] = array[1]

    protein_sequence = open(os.path.join(tmp_path, "protein.fasta"), "w")
    for value in longest_mRNA.values():
        protein_sequence.write(">" + value + "\n" + pro_seq[value] + "\n")
    protein_sequence.close()
    protein_file = os.path.join(tmp_path, "protein.fasta")
    return protein_file


def kaks_des(tmp_path, gff_file):
    longest_len = {}
    longest_mRNA = {}
    with open(os.path.join(tmp_path, "gene_trans_info.txt"), 'r') as f:
        lines = f.readlines()[1:]
        for line in lines:
            line = line.rstrip("\n")
            array = line.split("\t")
            if array[0] not in longest_len:
                longest_len[array[0]] = array[6]
                longest_mRNA[array[0]] = array[1]
            elif array[0] in longest_len and int(array[6]) > int(longest_len[array[0]]):
                longest_len[array[0]] = array[6]
                longest_mRNA[array[0]] = array[1]

    gene_site = {}
    with open(gff_file, 'r') as file:
        for l in file:
            l = l.rstrip("\n")
            array = l.split("\t")
            if len(array) > 7 and "gene" in array[2]:
                describe = re.split('[;,]', array[8])
                gene = re.split('[=]', describe[0])[-1]
                gene_site.setdefault(gene, []).append(array[0])  # chr
                gene_site.setdefault(gene, []).append(array[3])  # start
                gene_site.setdefault(gene, []).append(array[4])  # end

    mRNA_site = {}
    with open(os.path.join(tmp_path, "gene_trans_info.txt"), 'r') as file:
        ls = file.readlines()[1:]
        for l in ls:
            l = l.rstrip("\n")
            array = l.split("\t")
            if array[1] == longest_mRNA[array[0]]:
                mRNA_site.setdefault(array[0], []).append(array[1])  # mRNA
                mRNA_site.setdefault(array[0], []).append(array[2])  # chr
                mRNA_site.setdefault(array[0], []).append(array[4])  # start
                mRNA_site.setdefault(array[0], []).append(array[5])  # end
                mRNA_site.setdefault(array[0], []).append(array[7])  # cds_num
                mRNA_site.setdefault(array[0], []).append(array[8])  # cds_site

    des = {}
    strand = {}
    with open(os.path.join(tmp_path, "des_retro.txt"), 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            array = line.split("\t")
            strand[array[0]] = array[1]
            des[array[0]] = array[2]

    kaks_retro = {}
    with open(os.path.join(tmp_path, "kaks_result.txt"), 'r') as f:
        lines = f.readlines()[1:]
        for line in lines:
            line = line.rstrip("\n")
            array = line.split("\t")
            retro_name = array[0].split("&")[1]
            kaks_retro.setdefault(retro_name, []).append(array[2])  # Ka
            kaks_retro.setdefault(retro_name, []).append(array[3])  # Ks
            kaks_retro.setdefault(retro_name, []).append(array[4])  # Ka/Ks

    retro_seq = {}
    with open(os.path.join(tmp_path, "all_retro_dna.fa"), 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line.replace(">", "")
            else:
                retro_seq[name] = line.upper()

    retro_pro = {}
    with open(os.path.join(tmp_path, "all_retro_pro.fa"), 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line.replace(">", "")
            else:
                retro_pro[name] = line.upper()

    host_gene = {}
    retrogene_chimerical = {}
    OUT = open(os.path.join(tmp_path, "final.out"), 'w')
    OUT.write("Retrocopy_ID\tRetro_chr\tRetro_strand\tRetro_start\tRetro_end\tParental_gene_ID\tParent_chr\t" +
              "Parent_start\tParent_end\tPro_start\tPro_end\tIdentity\tCoverage\tLost_intron\tDescription\t" +
              "Ka\tKs\tKa/Ks\tHost_gene_ID\tRetro_sequence\tRetro_protein\n")
    with open(os.path.join(tmp_path, "best_parent.out"), 'r') as f:
        lines = f.readlines()[1:]
        for line in lines:
            line = line.rstrip("\n")
            array = line.split("\t")
            chr = array[1]
            start = int(array[2])
            end = int(array[3])
            half = (end - start) // 2
            host_gene[array[0]] = "NA"
            retrogene_chimerical[array[0]] = des[array[0]]
            for key in gene_site:
                if chr == gene_site[key][0] and int(end - half) >= int(gene_site[key][1]) \
                        and int(end) <= int(gene_site[key][2]):
                    host_gene[array[0]] = key
                elif chr == gene_site[key][0] and int(start) >= int(gene_site[key][1]) \
                        and int(start + half) <= int(gene_site[key][2]):
                    host_gene[array[0]] = key
                elif chr == gene_site[key][0] and int(start) <= int(gene_site[key][1]) \
                        and int(end) >= int(gene_site[key][2]):
                    host_gene[array[0]] = key
                    retrogene_chimerical[array[0]] = "retrogene"

            if host_gene[array[0]] in mRNA_site.keys() and retrogene_chimerical[array[0]] != "retrogene" :
                if int(mRNA_site[host_gene[array[0]]][4]) == 1:
                    if int(start + 10) >= int(mRNA_site[host_gene[array[0]]][2]) \
                            and int(end - 10) <= int(mRNA_site[host_gene[array[0]]][3]):
                        retrogene_chimerical[array[0]] = "retrogene"
                    else:
                        retrogene_chimerical[array[0]] = "chimerical"
                elif int(mRNA_site[host_gene[array[0]]][4]) > 1:
                    retrogene_chimerical[array[0]] = "chimerical"
                    cds_site = re.split('[-,]', mRNA_site[host_gene[array[0]]][5])
                    for site_i in range(1, len(cds_site)//2):
                        intron_start = cds_site[(2*site_i-1)]
                        intron_end = cds_site[2*site_i]
                        if int(start + 50) >= int(intron_start) and int(end -50) <= int(intron_end):
                            retrogene_chimerical[array[0]] = des[array[0]]
            if host_gene[array[0]] in longest_mRNA:
                host_gene[array[0]] = mRNA_site[host_gene[array[0]]][0]
            OUT.write(
                str(array[0]) + "\t" + str(array[1]) + "\t" + str(strand[array[0]]) + "\t" + str(array[2]) + "\t" +
                str(array[3]) + "\t" + str(array[4]) + "\t" + str(array[5]) + "\t" + str(array[6]) + "\t" +
                str(array[7]) + "\t" + str(array[8]) + "\t" + str(array[9]) + "\t" + str(array[10]) + "\t" +
                str(array[11]) + "\t" + str(array[12]) + "\t" + str(retrogene_chimerical[array[0]]) + "\t" +
                str(kaks_retro[array[0]][0]) + "\t" +
                str(kaks_retro[array[0]][1]) + "\t" + str(kaks_retro[array[0]][2]) + "\t" + host_gene[array[0]] + "\t"
                + retro_seq[array[0]] + "\t" + retro_pro[array[0]] + "\n")
    OUT.close()
    os.remove(os.path.join(tmp_path, "best_parent.out"))
    os.remove(os.path.join(tmp_path, "des_retro.txt"))
    os.remove(os.path.join(tmp_path, "kaks_result.txt"))
    os.remove(os.path.join(tmp_path, "all_retro_dna.fa"))
    os.remove(os.path.join(tmp_path, "all_retro_pro.fa"))


# generate gtf file of parental genes and retrocopies
def gtf_generate(tmp_path):
    parent_name = []
    retro_info = {}
    with open(os.path.join(tmp_path, "final.out"), 'r') as f:
        lines = f.readlines()[1:]
        for line in lines:
            line = line.rstrip("\n")
            array = line.split("\t")
            parent_name.append(array[5])
            retro_info.setdefault(array[0], []).append(array[1])  # chr
            retro_info.setdefault(array[0], []).append(array[2])  # strand
            retro_info.setdefault(array[0], []).append(array[3])  # start
            retro_info.setdefault(array[0], []).append(array[4])  # end

    OUT = open(os.path.join(tmp_path, "retro.gtf"), 'w')
    with open(os.path.join(tmp_path, "genome.gtf"), 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            array = line.split("\t")
            des = array[8].split(";")
            for j in range(len(des)):
                if "transcript_id" in des[j]:
                    mRNA = des[j].replace("transcript_id \"", "")
                    mRNA = mRNA.replace("\"", "")
            if mRNA in parent_name:
                OUT.write(line + "\n")

    for key in retro_info:
        OUT.write(str(retro_info[key][0]) + "\tRetrocopy\texon\t" + str(retro_info[key][2]) + "\t" +
                  str(retro_info[key][3]) + "\t.\t" + str(retro_info[key][1]) + "\t.\ttranscript_id \"" + str(key) +
                  "\"; gene_id \"" + str(key) + "\"; gene_name \"" + str(key) + "\";\n")
        OUT.write(str(retro_info[key][0]) + "\tRetrocopy\tCDS\t" + str(retro_info[key][2]) + "\t" +
                  str(retro_info[key][3]) + "\t.\t" + str(retro_info[key][1]) + "\t0\ttranscript_id \"" + str(key) +
                  "\"; gene_id \"" + str(key) + "\"; gene_name \"" + str(key) + "\";\n")

    OUT.close()


def START(args):
    genome_file = args['genome_file']
    gff_file = args['gff_file']
    working_dir = args['working_dir']

    gffread = args['gffread']
    lastdb = args['lastdb']
    lastal = args['lastal']
    bedtools = args['bedtools']

    identity = args['identity']
    coverage_rate = args['coverage_rate']
    coverage_len = args['coverage_len']
    intron_loss_num = args['intron_loss_num']
    gaplen = args['gaplen']
    deduplication_identity_cutoff = args['deduplication_identity_cutoff']
    deduplication_coverage_cutoff = args['deduplication_coverage_cutoff']
    deduplication_similar_cutoff = args['deduplication_similar_cutoff']
    thread = args['thread']
    kaksmethod = args['kaksmethod']

    seqkit = args['seqkit']
    diamond = args['diamond']
    clustalw2 = args['clustalw2']

    RNASeq_path = args['RNASeq_path']
    hisat2 = args['hisat2']
    stringtie = args['stringtie']
    samtools = args['samtools']
    hisat2build = args['hisat2build']

    bin_path = os.path.abspath(sys.argv[0]).replace("retroscan.py","")
    # 1. check input data correctness and set temp_path and logfile
    genome_file, gff_file, working_dir = check_input_data(genome_file, gff_file, working_dir)
    tmp_path = os.path.join(working_dir, "retroscan_results")
    logfile = os.path.join(tmp_path, "RetroScan.log")
    log = Logger(logfile, level='debug')

    # 2. Preparetion before alignment.
    log.logger.info("Preparetion before alignment.")
    log.logger.info("Use gffread to extract protein sequences and cds sequences.")
    gffread_cmd = "time bash " + os.path.join(bin_path, "pep_cds_generation.sh") + " " + tmp_path + " " + gffread + " " \
                  + gff_file + " " + genome_file + " " + logfile
    callCMD(gffread_cmd)
    if os.path.exists(os.path.join(tmp_path, "pep.fasta")) and os.path.getsize(os.path.join(tmp_path, "pep.fasta")) > 0\
            and os.path.exists(os.path.join(tmp_path, "cds.fasta")) and os.path.getsize(os.path.join(tmp_path, "cds.fasta")) > 0 :
        log.logger.info("Successfully extract protein sequences and cds sequences.")
    else:
        log.logger.error("Failed to extract protein sequences and cds sequences, please check gffread software or input files.")
        sys.exit(0)
    log.logger.info("Extract longest protein sequence for each gene which has three or more exons.")
    protein_path = os.path.join(tmp_path, "pep.fasta")
    protein_file = gene_trans_info(tmp_path, protein_path, gff_file)
    if os.path.exists(protein_file) and os.path.exists(os.path.join(tmp_path, "gene_trans_info.txt")) \
            and os.path.getsize(protein_file) > 0 and os.path.getsize(os.path.join(tmp_path, "gene_trans_info.txt")) > 0:
        log.logger.info("Preparetion has done.")
    else:
        log.logger.error("Preparetion has something wrong, please check input files, especially whether the gff3 format is standard.")
        sys.exit(0)

    # 3.algment and find best parent gene
    # 4.delete pairs which one parental gene vs a lot of retrocopies
    log.logger.info("Alignment.")
    alignment_cmd = "time bash " + os.path.join(bin_path, "alignment.sh") + " " + tmp_path + " " + protein_file + " " \
                    + genome_file + " " + lastdb + " " + lastal + " " + bedtools + " " + str(identity) + " " + \
                    str(coverage_rate) + " " + str(coverage_len) + " " + str(intron_loss_num)  + " " + str(gaplen) + " " \
                    + str(deduplication_identity_cutoff) + " " + str(deduplication_coverage_cutoff) + " " \
                    + str(deduplication_similar_cutoff) + " " + str(thread) + " " + logfile
    callCMD(alignment_cmd)
    if os.path.exists(os.path.join(tmp_path, "best_parent.out")) and os.path.getsize(os.path.join(tmp_path, "best_parent.out")) > 0 :
        log.logger.info("Successfully alignment and find best parental genes.")
    else:
        log.logger.error("Failed to align, please check LAST software.")
        sys.exit(0)

    # 5.calculate kaks and identity intact retrocopy, chimerical retrocopy, pseudogene
    log.logger.info("Calculate kaks and identity intact retrocopy, chimerical retrocopy, pseudogene.")
    kaks_intact_cmd = "time bash " + os.path.join(bin_path, "kaks_intact.sh") + " " + tmp_path + " " + protein_file \
                      + " " + genome_file + " " + bedtools + " " + kaksmethod + " " + " " + seqkit + " " + diamond + " " \
                      + clustalw2 + " " + str(thread) + " " + logfile
    callCMD(kaks_intact_cmd)
    if os.path.exists(os.path.join(tmp_path, "des_retro.txt")) and os.path.getsize(os.path.join(tmp_path, "des_retro.txt")) > 0 \
            and os.path.exists(os.path.join(tmp_path, "kaks_result.txt")) and os.path.getsize(os.path.join(tmp_path, "kaks_result.txt")) > 0 \
            and os.path.exists(os.path.join(tmp_path, "all_retro_dna.fa")) and os.path.getsize(os.path.join(tmp_path, "all_retro_dna.fa")) > 0:
        log.logger.info("Successfully calculate KaKs.")
        kaks_des(tmp_path, gff_file)
    else:
        log.logger.error("Calculate KaKs failed.")
        sys.exit(0)
    if os.path.exists(os.path.join(tmp_path, "final.out")) and os.path.getsize(os.path.join(tmp_path, "final.out")) > 0:
        log.logger.info("Successfully got the retrocopies results.")
    else:
        log.logger.error("Failed to get the retrocopies results.")
        sys.exit(0)

    # 6.Reprot
    final_out = pd.read_table(os.path.join(tmp_path, "final.out"), sep='\t', delimiter=None, header=0)

    with PdfPages(os.path.join(tmp_path, "Result_report.pdf")) as pdf:
        # The number distribution of retrocopies owned by each parent gene
        df = pd.DataFrame(final_out['Parental_gene_ID'].value_counts())
        num_df = pd.DataFrame(df['Parental_gene_ID'].value_counts())
        X = num_df._stat_axis.values.tolist()
        Y = num_df["Parental_gene_ID"]
        plt.bar(X, Y, 0.4)
        plt.xlabel("Number of parental genes")
        plt.ylabel("Number of retrocopies")
        plt.title("The number distribution of retrocopies owned by each parent gene")
        pdf.savefig()
        plt.close()

        # Distribution of retrocopy length
        final_out["length"] = final_out["Retro_end"] - final_out["Retro_start"] + 1
        bins = np.linspace(min(final_out["length"]), max(final_out["length"]), 100)
        plt.hist(final_out["length"], bins)
        plt.xlabel('Length of retrocopies')
        plt.ylabel('Number of retrocopies')
        plt.title('Distribution of retrocopy length')
        pdf.savefig()
        plt.close()

        # Percentage of identity
        labels = ["50% <= Identity < 60%", "60% <= Identity < 70%", "70% <= Identity < 80%", "80% <= Identity < 90%",
                  "90% <= Identity <= 100%"]
        x1 = final_out[(final_out['Identity'] <= 100) & (final_out['Identity'] >= 90)].iloc[:, 0].size
        x2 = final_out[(final_out['Identity'] < 90) & (final_out['Identity'] >= 80)].iloc[:, 0].size
        x3 = final_out[(final_out['Identity'] < 80) & (final_out['Identity'] >= 70)].iloc[:, 0].size
        x4 = final_out[(final_out['Identity'] < 70) & (final_out['Identity'] >= 60)].iloc[:, 0].size
        x5 = final_out[(final_out['Identity'] < 60) & (final_out['Identity'] >= 50)].iloc[:, 0].size
        X = [x1, x2, x3, x4, x5]
        plt.pie(X, labels=labels, colors=('yellow', 'lightskyblue', 'lightgreen', 'pink', 'sandybrown'),
                autopct='%1.1f%%')
        plt.title("Percentage of identity")
        pdf.savefig()
        plt.close()

        # Percentage of Coverage
        labels = ["50% <= Coverage < 60%", "60% <= Coverage < 70%", "70% <= Coverage < 80%", "80% <= Coverage < 90%",
                  "90% <= Coverage <= 100%"]
        x1 = final_out[(final_out['Coverage'] <= 100) & (final_out['Coverage'] >= 90)].iloc[:, 0].size
        x2 = final_out[(final_out['Coverage'] < 90) & (final_out['Coverage'] >= 80)].iloc[:, 0].size
        x3 = final_out[(final_out['Coverage'] < 80) & (final_out['Coverage'] >= 70)].iloc[:, 0].size
        x4 = final_out[(final_out['Coverage'] < 70) & (final_out['Coverage'] >= 60)].iloc[:, 0].size
        x5 = final_out[(final_out['Coverage'] < 60) & (final_out['Coverage'] >= 50)].iloc[:, 0].size
        X = [x1, x2, x3, x4, x5]
        plt.pie(X, labels=labels, colors=('yellow', 'lightskyblue', 'lightgreen', 'pink', 'sandybrown'),
                autopct='%1.1f%%')
        plt.title("Percentage of coverage")
        pdf.savefig()
        plt.close()

        # Percentage of category
        labels = ["retrogene", "intact", "chimerical", "pseudogene"]
        x1 = final_out[final_out['Description'] == "retrogene"].iloc[:, 0].size
        x2 = final_out[final_out['Description'] == "intact"].iloc[:, 0].size
        x3 = final_out[final_out['Description'] == "chimerical"].iloc[:, 0].size
        x4 = final_out[final_out['Description'] == "pseudogene"].iloc[:, 0].size
        X = [x1, x2, x3, x4]
        plt.pie(X, labels=labels, colors=('yellow', 'lightskyblue', 'lightgreen', 'pink'), autopct='%1.1f%%')
        plt.title("Percentage of category")
        pdf.savefig()
        plt.close()

        # Distribution of Ks
        ks_len = final_out[final_out['Ks'] != "NA"]
        bins = np.linspace(min(ks_len['Ks']), max(ks_len['Ks']), 100)
        plt.hist(ks_len['Ks'], bins)
        plt.xlabel('Ks')
        plt.ylabel('Number of retrocopies')
        plt.title('Distribution of Ks')
        pdf.savefig()
        plt.close()

        # Distribution of Ka/Ks
        kaks_len = final_out[final_out['Ka/Ks'] != "NA"]
        bins = np.linspace(min(kaks_len['Ka/Ks']), max(kaks_len['Ka/Ks']), 100)
        plt.hist(kaks_len['Ka/Ks'], bins)
        plt.xlabel('Ka/Ks')
        plt.ylabel('Number of retrocopies')
        plt.title('Distribution of Ka/Ks')
        pdf.savefig()
        plt.close()

    # 7.RNA-Seq
    if RNASeq_path == "NA" or RNASeq_path == "":
        log.logger.info("No RNA-Seq data.")
        log.logger.info("Done.")
    elif os.path.exists(os.path.abspath(RNASeq_path)):
        log.logger.info("RNA-Seq")
        gtf_cmd = gffread + " " + gff_file + " -T -o " + os.path.join(tmp_path, "genome.gtf")
        callCMD(gtf_cmd)
        gtf_generate(tmp_path)
        gtf_file = os.path.join(tmp_path, "retro.gtf")
        RNASeq_cmd = "time bash " + os.path.join(bin_path, "RNASeq.sh") + " " + os.path.abspath(RNASeq_path) + " " + \
                     genome_file + " " + gtf_file + " " + str(thread) + " " + hisat2 + " " + stringtie + " " + \
                     samtools + " " + hisat2build + " " + tmp_path + " " + logfile
        callCMD(RNASeq_cmd)
        os.remove(os.path.join(tmp_path, "genome.gtf"))
        os.remove(os.path.join(tmp_path, "retro.gtf"))
        if os.path.exists(os.path.join(tmp_path, "all_samples.counts.txt")) and os.path.getsize(os.path.join(tmp_path, "all_samples.counts.txt")) > 0:
            log.logger.info("Successfully calculate expression values.")
            log.logger.info("Done.")
        else:
            log.logger.error("Calculate expression values failed.")
            sys.exit(0)


def main():
    argv = sys.argv
    parser = argparse.ArgumentParser()

    parser.add_argument('genome_file', metavar='genome_file', type=str,
                        help='- Genome sequences fasta file, ending with fasta, fa, fas, faa, fna, tar.gz or zip')
    parser.add_argument('gff_file', metavar='gff_file', type=str,
                        help='- Gff file, ending with gff, gff3, tar.gz or zip')
    parser.add_argument('working_dir', metavar='Output_dir', type=str,
                        help='- Path to the directory where all intermediate and final results will be stored')

    parser.add_argument('--identity', metavar='int', type=int, nargs='?', default=50,
                        help='- Expectation identity threshold(percentage) for saving hits [50] ')
    parser.add_argument('--coverage_rate', metavar='int', type=int, nargs='?', default=50,
                        help='- Expectation coverage rate(percentage) threshold for saving hits [50] ')
    parser.add_argument('--coverage_len', metavar='int', type=int, nargs='?', default=50,
                        help='- Expectation length(aa) of coverage for saving hits [50] ')
    parser.add_argument('--intron_loss_num', metavar='int', type=int, nargs='?', default=2,
                        help='- Expectation the number of intron-lost for saving hits [2] ')
    parser.add_argument('--gaplen', metavar='int', type=int, nargs='?', default=40,
                        help='- Expectation length(bp) of gap for merging hits [40] ')
    parser.add_argument('--deduplication_identity_cutoff', metavar='int', type=int, nargs='?', default=80,
                        help='- The identity cutoff (percentage) of retrocopies duplication sequences aligned to the genome [80] ')
    parser.add_argument('--deduplication_coverage_cutoff', metavar='int', type=int, nargs='?', default=80,
                        help='- The coverage cutoff (percentage) of retrocopies duplication sequences aligned to the genome [80] ')
    parser.add_argument('--deduplication_similar_cutoff', metavar='int', type=int, nargs='?', default=10,
                        help='- The cutoff number of retrocopies duplication sequences aligned to the genome [10] ')
    parser.add_argument('--kaksmethod', type=str, nargs='?', default='NG', choices=['NG', 'LWL', 'LPB', 'MLWL', 'MLPB', 'GY', 'YN', 'MYN', 'MS', 'MA', 'GNG', 'GLWL', 'GLPB', 'GMLWL', 'GMLPB', 'GYN', 'GMYN'],
                        help="- Methods for estimating Ka and Ks and theirs references, choices=['NG', 'LWL', 'LPB', 'MLWL', 'MLPB', 'GY', 'YN', 'MYN', 'MS', 'MA', 'GNG', 'GLWL', 'GLPB', 'GMLWL', 'GMLPB', 'GYN', 'GMYN']")
    parser.add_argument('--thread', metavar='int', type=int, nargs='?', default=1,
                        help='- Number of threads (CPUs) to use [1] ')

    parser.add_argument('--gffread', type=str, nargs='?', default="gffread",
                        help='- Path to the gffread, https://github.com/gpertea/gffread (Default is "gffread").')
    parser.add_argument('--lastdb', type=str, nargs='?', default="lastdb",
                        help='- Path to the lastdb, http://last.cbrc.jp/ (Default is "lastdb").')
    parser.add_argument('--lastal', type=str, nargs='?', default="lastal",
                        help='- Path to the lastal, http://last.cbrc.jp/ (Default is "lastal").')
    parser.add_argument('--bedtools', type=str, nargs='?', default="bedtools",
                        help='- Path to the bedtools, https://bedtools.readthedocs.io/en/latest/index.html (Default is "bedtools").')

    parser.add_argument('--seqkit', type=str, nargs='?', default="seqkit",
                        help='- Path to the seqkit, https://bioinf.shenwei.me/seqkit/ (Default is "seqkit").')
    parser.add_argument('--diamond', type=str, nargs='?', default="diamond",
                        help='- Path to the diamond, http://github.com/bbuchfink/diamond (Default is "diamond").')
    parser.add_argument('--clustalw2', type=str, nargs='?', default="clustalw2",
                        help='- Path to the clustalw2, http://www.clustal.org/clustal2/ (Default is "clustalw2").')
 

    parser.add_argument('--RNASeq_path', type=str, nargs='?', default="",
                        help='- Path to the RNASeq file.')
    parser.add_argument('--hisat2', type=str, nargs='?', default="hisat2",
                        help='- Path to the hisat2, https://daehwankimlab.github.io/hisat2/ (Default is "hisat2").')
    parser.add_argument('--hisat2build', type=str, nargs='?', default="hisat2-build",
                        help='- Path to the hisat2-build, https://daehwankimlab.github.io/hisat2/ (Default is "hisat2-build").')
    parser.add_argument('--stringtie', type=str, nargs='?', default="stringtie ",
                        help='- Path to the stringtie, http://ccb.jhu.edu/software/stringtie/ (Default is "stringtie ").')
    parser.add_argument('--samtools', type=str, nargs='?', default="samtools",
                        help='- Path to the samtools, http://www.htslib.org/ (Default is "samtools").')


    parser.add_argument('--version', action='version', version='RetroScan version 1.0.0')


    args = vars(parser.parse_args())

    START(args)


main()
