#!/bin/bash
tmp_path=$1
gffread=$2
gff_file=$3
genome_file=$4
LOGFILE=$5

pipeline() {
if [[ ! -x $gffread ]]; then
	gffread=`which gffread`
	if [[ ! -x $gffread ]]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: gffread program not found."
	exit 1
	fi
fi

if [ ! -d $tmp_path ]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The result direction is wrong."
	exit 1
fi

if [ ! -e $genome_file ]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The genome file does not exist."
	exit 1
fi

if [ ! -e $gff_file ]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The gff file does not exist."
	exit 1
fi

if [ ! -e $LOGFILE ]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: The log file does not exist."
	exit 1
fi


echo [`date +"%Y-%m-%d %H:%M:%S"`] " Begin to generate cds.fasta and pep.fasta."
$gffread $gff_file -g $genome_file -x $tmp_path/cds.fasta -y $tmp_path/tmp_pep.fasta
sed -i 's/\.$//g' $tmp_path/tmp_pep.fasta
sed -i '/^\s*$/d' $tmp_path/tmp_pep.fasta
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\t":$0 }' $tmp_path/tmp_pep.fasta | awk -F'\t' '{if($2!~"\\."){print $1"\n"$2}}' > $tmp_path/pep.fasta
rm $tmp_path/tmp*
echo [`date +"%Y-%m-%d %H:%M:%S"`] " End to generate cds.fasta and pep.fasta."
}

pipeline 2>&1 | tee -a $LOGFILE
