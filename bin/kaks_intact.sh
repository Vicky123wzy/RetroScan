#!/bin/bash
tmp_path=$1
protein_file=$2
genome_file=$3
bedtools=$4
kaksmethod=$5

seqkit=$6
diamond=$7
clustalw2=$8

thread=$9
LOGFILE=${10}

kaks_intact () {
bin_path=`dirname $0`

if [[ ! -x $bedtools ]]; then
	bedtools=`which bedtools`
	if [[ ! -x $bedtools ]]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: bedtools program does not found."
	exit 1
	fi
fi

if [[ ! -x $seqkit ]]; then
	seqkit=`which seqkit`
	if [[ ! -x $seqkit ]]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: seqkit program does not found."
	exit 1
	fi
fi

if [[ ! -x $diamond ]]; then
	diamond=`which diamond`
	if [[ ! -x $diamond ]]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: diamond program does not found."
	exit 1
	fi
fi

if [[ ! -x $clustalw2 ]]; then
	clustalw2=`which clustalw`
	if [[ ! -x $clustalw2 ]]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: clustalw2 program does not found."
	exit 1
	fi
fi


if [ ! -d $tmp_path ]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The file path does exist."
	exit 1
fi


python $bin_path/kaks_thread.py $tmp_path/best_parent.out $tmp_path $bin_path/kaks_single.sh $genome_file $protein_file $bedtools $kaksmethod $seqkit $diamond $clustalw2 $thread

cat $tmp_path/kaks_result/*/retro_dna.fasta > $tmp_path/all_retro_dna.fa
cat $tmp_path/kaks_result/*/retro_pro.fasta > $tmp_path/all_retro_pro.fa
cat $tmp_path/kaks_result/*/des_retro.txt > $tmp_path/des_retro.txt
cat $tmp_path/kaks_result/*/kaks_result.txt > $tmp_path/kaks_result.txt

sed -i '1i\Sequence\tMethod\tKa\tKs\tKa/Ks\tP-Value(Fisher)\tLength\tS-Sites\tN-Sites\tFold-Sites(0:2:4)\tSubstitutions\tS-Substitutions\tN-Substitutions\tFold-S-Substitutions(0:2:4)\tFold-N-Substitutions(0:2:4)\tDivergence-Time\tSubstitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)\tGC(1:2:3)\tML-Score\tAICc\tAkaike-Weight\tModel' $tmp_path/kaks_result.txt

rm -rf $tmp_path/kaks_result
echo [`date +"%Y-%m-%d %H:%M:%S"`] " Finish the kaks calculator."
}
kaks_intact  2>&1 | tee -a $LOGFILE

