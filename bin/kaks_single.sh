#!/bin/bash
tmp_path=$1
genome_file=$2
protein_file=$3
lossintronfile=$4
cdsfile=$5
bedtools=$6
kaksmethod=$7

seqkit=$8
diamond=$9
clustalw2=${10}

kaks_intact () {
bin_path=`dirname $0`
AXTConvertor=$bin_path/AXTConvertor
KaKs_Calculator=$bin_path/KaKs_Calculator

retro=$(cat $lossintronfile | awk '{print $1}')
cat $lossintronfile | awk '{print $2"\t"$3-1"\t"$4}' > $tmp_path/temp_retro_dna.bed
$bedtools getfasta -fi $genome_file -bed $tmp_path/temp_retro_dna.bed > $tmp_path/retro_dna.fasta
sed -i '1c >'$retro'' $tmp_path/retro_dna.fasta

$seqkit translate -F -f 6 $tmp_path/retro_dna.fasta | sed -r '/^6/d' > $tmp_path/temp_retro_pro.fasta
mRNA=$(cat $lossintronfile | awk '{print ">"$5}')
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\t":$0 }' $protein_file | awk '{if($1=="'"$mRNA"'"){print $1"\n"$NF}}' > $tmp_path/temp_parent_pro.fasta
sed -i '1c '$mRNA'' $tmp_path/temp_parent_pro.fasta
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\t":$0 }' $cdsfile | awk '{if($1=="'"$mRNA"'"){print $1"\n"$NF}}' > $tmp_path/temp_parent_dna.fasta
sed -i '1c '$mRNA'' $tmp_path/temp_parent_dna.fasta
$diamond makedb -p 1 --in $tmp_path/temp_parent_pro.fasta -d $tmp_path/temp_parent
$diamond blastp -d $tmp_path/temp_parent -q $tmp_path/temp_retro_pro.fasta --more-sensitive --max-target-seqs 1000 -f 6 -p 1 -o $tmp_path/temp_blast.out

python $bin_path/remove_stop_codon.py $tmp_path/temp_blast.out $tmp_path/temp_retro_pro.fasta $tmp_path/retro_dna.fasta $tmp_path/temp_retro_dna1.fasta $tmp_path/temp_retro_pro1.fasta $tmp_path/des_retro.txt $tmp_path/retro_pro.fasta
sed -i '1c >'$retro'' $tmp_path/retro_pro.fasta

cat $tmp_path/temp_parent_dna.fasta $tmp_path/temp_retro_dna1.fasta > $tmp_path/temp_dna.fasta
cat $tmp_path/temp_parent_pro.fasta $tmp_path/temp_retro_pro1.fasta > $tmp_path/temp_pro.fasta
$clustalw2 -INFILE=$tmp_path/temp_pro.fasta -OUTPUT=CLUSTAL -OUTFILE=$tmp_path/temp_pro.aln
perl $bin_path/pal2nal.pl $tmp_path/temp_pro.aln $tmp_path/temp_dna.fasta -output clustal -nogap > $tmp_path/temp_codon
echo $AXTConvertor
echo $KaKs_Calculator
$AXTConvertor $tmp_path/temp_codon $tmp_path/temp.axt
$KaKs_Calculator -i $tmp_path/temp.axt -o $tmp_path/temp_kaks.txt -m $kaksmethod
awk 'NR==2{print}' $tmp_path/temp_kaks.txt > $tmp_path/kaks_result.txt

rm $tmp_path/temp* 

}
kaks_intact  2>&1 | tee -a $LOGFILE
