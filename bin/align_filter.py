import sys
import re
import os

USAGE = "\nusage: python %s inputfile outputfile gene_trans_info identity coverage_rate coverage_len intron_loss_num\n" % sys.argv[0]

if len(sys.argv) != 9:
    print(USAGE)
    sys.exit()

inputfile = sys.argv[1]
outputfile = sys.argv[2]
gene_trans_info = sys.argv[3]
cluster_file = sys.argv[4]

identity = sys.argv[5]
coverage_rate = sys.argv[6]
coverage_len = sys.argv[7]
intron_loss_num = sys.argv[8]

pep_site = {}
mRNA_info = {}
with open(gene_trans_info, 'r') as f:
    for line in f:
        line = line.rstrip("\n")
        array = line.split("\t")
        mRNA_info.setdefault(array[1], []).append(array[2])  # chr
        mRNA_info.setdefault(array[1], []).append(array[4])  # start
        mRNA_info.setdefault(array[1], []).append(array[5])  # end
        mRNA_info.setdefault(array[1], []).append(array[6])  # protein_length
        pep_site[array[1]] = array[10]

cluster_num = 0
best_seq = {}
with open(cluster_file, 'r') as file:
    for line in file:
        line = line.rstrip("\n")
        array = line.split("\t")
        cluster_num = cluster_num + 1
        best_seq.setdefault(cluster_num, []).append(array[0]) #chr
        best_seq.setdefault(cluster_num, []).append(array[1]) #start
        best_seq.setdefault(cluster_num, []).append(array[2]) #end
        best_seq.setdefault(cluster_num, []).append("0") #best_score
        best_seq.setdefault(cluster_num, []).append("NA") #line

with open(inputfile, 'r') as file:
    for line in file:
        line = line.rstrip("\n")
        array = line.split("\t")
        if (
                (len(array)==12 and array[0] != mRNA_info[array[1]][0])
                or (len(array)==12 and array[0] == mRNA_info[array[1]][0] and int(array[6]) <= int(mRNA_info[array[1]][1])
                    and int(array[7]) <= int(mRNA_info[array[1]][1]))
                or (len(array)==12 and array[0] == mRNA_info[array[1]][0] and int(array[6]) >= int(mRNA_info[array[1]][2])
                    and int(array[7]) >= int(mRNA_info[array[1]][2]))
        ):
            parent_pro = array[1]
            identity_score = array[2]
            hit_cover_len = abs(int(array[9]) - int(array[8])) + 1
            hit_cover_rate = hit_cover_len / int(mRNA_info[parent_pro][3]) * 100
            intron_site = pep_site[parent_pro].split(",")
            selected = [x for x in intron_site if
                        int(x) in range(min(int(array[8]), int(array[9]))+10, max(int(array[8]), int(array[9]))-10)]
            intron_loss = len(selected)
            if float(identity_score) >= float(identity) and int(hit_cover_len) >= int(coverage_len) \
                    and float(hit_cover_rate) >= float(coverage_rate) and int(intron_loss) >= int(intron_loss_num):
                for key in best_seq:
                    if array[0] == best_seq[key][0] and min(int(array[6]),int(array[7])) >= int(best_seq[key][1])\
                            and max(int(array[6]),int(array[7])) <= int(best_seq[key][2]) and float(array[11]) > float(best_seq[key][3]):
                        best_seq[key][3] = array[11]
                        array.append(hit_cover_rate)
                        array.append(intron_loss)
                        best_seq[key][4] = array


output_file = open(outputfile, "w")
output_file.write("Retrocopy_ID\tRetro_chr\tRetro_start\tRetro_end\tParental_gene_ID\tParent_chr\tParent_start\t" +
                      "Parent_end\tPro_start\tPro_end\tIdentity\tCoverage\tLost_intron\n")
i=0
for key in best_seq:
    if best_seq[key][4] != "NA":
        array = best_seq[key][4]
        i = i + 1
        output_file.write("Retrocopy" + str(i) + "\t" + str(array[0]) + "\t" + str(min(int(array[6]), int(array[7])))
                          + "\t" + str(max(int(array[6]), int(array[7]))) + "\t" + str(array[1]) + "\t" +
                          str(mRNA_info[array[1]][0]) + "\t" + str(mRNA_info[array[1]][1]) + "\t" +
                          str(mRNA_info[array[1]][2]) + "\t" + str(min(int(array[8]), int(array[9]))) + "\t" +
                          str(max(int(array[8]), int(array[9]))) + "\t" + str(array[2]) + "\t" +
                          str(round(array[12], 1)) + "\t" + str(array[13]) + "\n")
output_file.close()
