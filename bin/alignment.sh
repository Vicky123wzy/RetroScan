#!/bin/bash
tmp_path=$1
protein_file=$2
genome_file=$3

#alignment_tool=$(grep "^alignment_tool" $config_file | sed 's/ //g' - | awk -F'=' '{print $2}' -)
lastdb=$4
lastal=$5
bedtools=$6

identity=$7
coverage_rate=$8
coverage_len=$9
intron_loss_num=${10}
gaplen=${11}
deduplication_identity_cutoff=${12}
deduplication_coverage_cutoff=${13} 
deduplication_similar_cutoff=${14}

thread=${15}
LOGFILE=${16}

pipeline() {
bin_path=`dirname $0`

if [ ! -d $tmp_path ]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The work path does not exist."
	exit 1
fi
if [ ! -d $bin_path ]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The bin path does not exist."
	exit 1
fi

if [[ ! -x $lastdb ]]; then
	lastdb=`which lastdb`
	if [[ ! -x $lastdb ]]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: LAST program does not found."
	exit 1
	fi
fi

if [[ ! -x $lastal ]]; then
	lastal=`which lastal`
	if [[ ! -x $lastal ]]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: LAST program does not found."
	exit 1
	fi
fi

if [[ ! -x $bedtools ]]; then
	bedtools=`which bedtools`
	if [[ ! -x $bedtools ]]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: bedtools program does not found."
	exit 1
	fi
fi


if [ ! -e $genome_file ]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The genome file does not exist."
	exit 1
fi
if [ ! -e $protein_file ]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The protein file does not exist."
	exit 1
fi
if [ ! -e $LOGFILE ]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The log file does not exist."
	exit 1
fi

if [[ $(echo "$identity > 0" | bc) = 1 && $(echo "$identity < 100" | bc) = 1 ]];then 
	echo " "
else
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The identity setting is wrong, identity must be between 0 and 100."
	exit 1
fi

if [[ $(echo "$coverage_rate > 0" | bc) = 1 && $(echo "$coverage_rate < 100" | bc) = 1 ]];then 
	echo " "
else
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The coverage rate setting is wrong, coverage rate must be between 0 and 100."
	exit 1
fi

if [[ $(echo $coverage_len/1|bc) != "$coverage_len" ]]; then 
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The coverage length setting is wrong, coverage length must be an integer."
	exit 1
fi

if [[ $(echo $gaplen/1|bc) != "$gaplen" ]]; then 
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The gap length setting is wrong, gap length must be an integer."
	exit 1
fi

if [[ $(echo "$deduplication_identity_cutoff > 0" | bc) = 1 && $(echo "deduplication_identity_cutoff < 100" | bc) = 1 ]];then 
	echo " "
else
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The deduplication_identity_cutoff setting is wrong, identity must be between 0 and 100."
	exit 1
fi

if [[ $(echo "$deduplication_coverage_cutoff > 0" | bc) = 1 && $(echo "$deduplication_coverage_cutoff < 100" | bc) = 1 ]];then 
	echo " "
else
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The deduplication_coverage_cutoff setting is wrong, coverage rate must be between 0 and 100."
	exit 1
fi

if [[ $(echo $deduplication_similar_cutoff/1|bc) != "$deduplication_similar_cutoff" ]]; then 
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The deduplication_similar_cutoff setting is wrong, it must be an integer."
	exit 1
fi

if [[ $(echo $thread/1|bc) != "$thread" ]]; then 
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The thread setting is wrong, thread must be an integer."
	exit 1
fi



echo [`date +"%Y-%m-%d %H:%M:%S"`] " Begin to align protein sequences with genome sequences."
$lastdb -P $thread -p -cR01 $tmp_path/temp_protdb $protein_file
$lastal -p BL62 -a 11 -b 2 -F 15 -x 20 -d 63 -z 20 -P $thread $tmp_path/temp_protdb $genome_file -f BlastTab > $tmp_path/last_result.txt
awk -F"\t" 'BEGIN{OFS="\t"} {if($7<=$8){print $1"\t"$7"\t"$8}else {print $1"\t"$8"\t"$7}}' $tmp_path/last_result.txt > $tmp_path/temp_last_align_filter.out
	
$bedtools sort -i $tmp_path/temp_last_align_filter.out | $bedtools merge -i - -d $gaplen > ${tmp_path}/temp_cluster_result.bed

python $bin_path/align_filter.py $tmp_path/last_result.txt $tmp_path/temp_best_parent.out ${tmp_path}/gene_trans_info.txt ${tmp_path}/temp_cluster_result.bed $identity $coverage_rate $coverage_len $intron_loss_num
	
awk 'NR != 1 {print $2"\t"$3"\t"$4}' $tmp_path/temp_best_parent.out > $tmp_path/temp_filter.bed
$bedtools getfasta -fi $genome_file -bed $tmp_path/temp_filter.bed > $tmp_path/temp_all_retro.fasta
$lastdb -uNEAR -R01  -P $thread $tmp_path/temp_genomedb $genome_file
$lastal $tmp_path/temp_genomedb $tmp_path/temp_all_retro.fasta -P $thread -f BlastTab > $tmp_path/temp_filter_last_result.txt
python $bin_path/deduplication.py $tmp_path/temp_filter_last_result.txt $tmp_path/temp_best_parent.out $tmp_path/best_parent.out $deduplication_identity_cutoff $deduplication_coverage_cutoff $deduplication_similar_cutoff


rm ${tmp_path}/temp*
rm ${tmp_path}/last_result.txt
rm ${tmp_path}/pep.fasta

echo [`date +"%Y-%m-%d %H:%M:%S"`] " End to align."

}

pipeline 2>&1 | tee -a $LOGFILE



