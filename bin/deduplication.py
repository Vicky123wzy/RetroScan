#!/usr/bin/env python
import sys
import re
import os

USAGE = "\nusage: python %s temp_filter_last_result alignment_result output_file score_cutoff similar_cutoff \n" % sys.argv[0]

if len(sys.argv) != 7:
    print(USAGE)
    sys.exit()

temp_filter_last_result = sys.argv[1]
alignment_result = sys.argv[2]
output_file = sys.argv[3]

identity_cutoff = sys.argv[4]
coverage_cutoff = sys.argv[5]
similar_cutoff = sys.argv[6]


match_num = {}
with open(temp_filter_last_result, 'r') as f:
    for line in f:
        line = line.rstrip("\n")
        array = line.split("\t")
        if len(array) == 12:
            query_length = re.split("[:-]",array[0])
            coverage = (abs(int(array[7]) - int(array[6])) + 1) / (int(query_length[2]) - int(query_length[1])) * 100
            if array[0] not in match_num and float(array[2]) >= float(identity_cutoff) \
                    and float(coverage) >= float(coverage_cutoff):
                match_num[array[0]] = 1
            if array[0] in match_num and float(array[2]) >= float(identity_cutoff) \
                    and float(coverage) >= float(coverage_cutoff):
                match_num[array[0]] = match_num[array[0]] + 1

OUT = open(output_file, "w")
OUT.write("Retrocopy_ID\tRetro_chr\tRetro_start\tRetro_end\tParental_gene_ID\tParent_chr\tParent_start\tParent_end\t" +
          "Pro_start\tPro_end\tIdentity\tCoverage\tLost_intron\n")
i = 0
with open(alignment_result, 'r') as f:
    lines = f.readlines()[1:]
    for line in lines:
        line = line.rstrip("\n")
        array = line.split("\t")
        last_query_id = array[1] + ":" + array[2] + "-" + array[3]
        if last_query_id in match_num and (match_num[last_query_id]) <= int(similar_cutoff):
            i = i + 1
            array[0] = "Retrocpy" + str(i)
            OUT.write("\t".join(str(i) for i in array) + "\n")
OUT.close()

