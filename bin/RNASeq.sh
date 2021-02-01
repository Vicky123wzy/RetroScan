#!/usr/bin/env bash
FASTQLOC=$1 
genome=$2
GTFFILE=$3 
NUMCPUS=$4 
HISAT2=$5 
STRINGTIE=$6
SAMTOOLS=$7 
hisat2build=$8 
WRKDIR=$9
LOGFILE=${10}
pipeline() {

if [[ ! -x $SAMTOOLS ]]; then
	SAMTOOLS=`which samtools`
	if [[ ! -x $SAMTOOLS ]]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: samtools program not found."
	exit 1
	fi
fi
if [[ ! -x $HISAT2 ]]; then
	HISAT2=`which hisat2`
	if [[ ! -x $HISAT2 ]]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: hisat2 program not found."
	exit 1
	fi
fi
if [[ ! -x $STRINGTIE ]]; then
	STRINGTIE=`which stringtie`
	if [[ ! -x $STRINGTIE ]]; then
	echo [`date +"%Y-%m-%d %H:%M:%S"`] "ERROR: stringtie program not found."
	exit 1
	fi
fi


ALIGNLOC=${WRKDIR}/hisat2
if [ ! -d $ALIGNLOC ]; then
   mkdir -p $ALIGNLOC
fi

BALLGOWNLOC=${WRKDIR}/ballgown
if [ ! -d $BALLGOWNLOC ]; then
   mkdir -p $BALLGOWNLOC
fi

# main script block

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START " 

$hisat2build $genome ${ALIGNLOC}/genomeindex

reads=(${FASTQLOC}/*)
reads=("${reads[@]##*/}") 

for ((i=0; i<=${#reads[@]}-1; i ++ )); do
    thisreads=("${reads[$i]%%.*}") 
    if [[ $thisreads =~ _1$ ]]; then
       sample="${thisreads%_*}"  
       reads1=${reads[$i]} 
       reads2=("${reads1/_1./_2.}") 
       stime=`date +"%Y-%m-%d %H:%M:%S"`
       echo "[$stime] Processing sample: $sample"
       echo [$stime] "   *Pair-end Alignment of reads to genome (HISAT2)"
       $HISAT2 -p $NUMCPUS --dta -x ${ALIGNLOC}/genomeindex \
        -1 ${FASTQLOC}/${reads1} \
        -2 ${FASTQLOC}/${reads2} \
        -S ${ALIGNLOC}/${sample}.sam 2>${ALIGNLOC}/${sample}.alnstats
       echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Alignments conversion (SAMTools)"
       $SAMTOOLS view -S -b ${ALIGNLOC}/${sample}.sam | \
       $SAMTOOLS sort -@ $NUMCPUS -o ${ALIGNLOC}/${sample}.bam -
   
       #echo "..removing intermediate files"
       rm ${ALIGNLOC}/${sample}.sam
       #rm ${TEMPLOC}/${sample}.unsorted.bam

       echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Assemble transcripts (StringTie)"
       $STRINGTIE -p $NUMCPUS -G ${GTFFILE} -o ${ALIGNLOC}/${sample}.gtf -l ${sample} ${ALIGNLOC}/${sample}.bam
    #   $HTSEQ -f bam -r name -s no -n $NUMCPUS -i transcript_id ${ALIGNLOC}/${sample}.bam $GTFFILE > ${ALIGNLOC}/${sample}.counts.txt
    elif [[ $thisreads =~ _2$ ]]; then
       echo " " 
    elif [[ $thisreads != _1$ ]]&&[[ $thisreads != _2$ ]]; then
       sample=$thisreads
       stime=`date +"%Y-%m-%d %H:%M:%S"`
       echo "[$stime] Processing sample: $sample"
       echo [$stime] "   *Single-end Alignment of reads to genome (HISAT2)"
       $HISAT2 -p $NUMCPUS --dta -x ${ALIGNLOC}/genomeindex \
        -U ${FASTQLOC}/${reads[$i]} \
        -S ${ALIGNLOC}/${sample}.sam 2>${ALIGNLOC}/${sample}.alnstats
       echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Alignments conversion (SAMTools)"
       $SAMTOOLS view -S -b ${ALIGNLOC}/${sample}.sam | \
       $SAMTOOLS sort -@ $NUMCPUS -o ${ALIGNLOC}/${sample}.bam -
   
       #echo "..removing intermediate files"
       rm ${ALIGNLOC}/${sample}.sam
       #rm ${TEMPLOC}/${sample}.unsorted.bam

       echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Assemble transcripts (StringTie)"
       $STRINGTIE -p $NUMCPUS -G ${GTFFILE} -o ${ALIGNLOC}/${sample}.gtf -l ${sample} ${ALIGNLOC}/${sample}.bam
     #  $HTSEQ -f bam -r name -s no -n $NUMCPUS -i transcript_id ${ALIGNLOC}/${sample}.bam $GTFFILE > ${ALIGNLOC}/${sample}.counts.txt
    fi
done

## merge transcript file
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Merge all transcripts (StringTie)"
ls -1 ${ALIGNLOC}/*.gtf > ${ALIGNLOC}/mergelist.txt

$STRINGTIE --merge -p $NUMCPUS -G  ${GTFFILE} \
    -o ${ALIGNLOC}/stringtie_merged.gtf ${ALIGNLOC}/mergelist.txt

## estimate transcript abundance
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Estimate abundance for each sample (StringTie)"
reads=(${ALIGNLOC}/*.bam)
reads=("${reads[@]##*/}")
for ((i=0; i<=${#reads[@]}-1; i ++ )); do
    sample="${reads[$i]%.*}"
    if [ ! -d ${BALLGOWNLOC}/${sample} ]; then
       mkdir -p ${BALLGOWNLOC}/${sample}
    fi
    $STRINGTIE -e -B -p $NUMCPUS -G ${ALIGNLOC}/stringtie_merged.gtf \
    -o ${BALLGOWNLOC}/${sample}/${sample}.gtf ${ALIGNLOC}/${sample}.bam
done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Generate the expression matrix (FPKM)"
for i in ${BALLGOWNLOC}/*; do
	samplename=("${i##*/}") 
	awk -F"\t" '{print $6"\t"$12}' $i/t_data.ctab > $i/fpkm.txt 
	sed -i '1c\gene\t'$samplename'' $i/fpkm.txt 
done

paste ${BALLGOWNLOC}/*/fpkm.txt | awk '{printf $1;for(i=2;i<=NF;i=i+2)printf "\t"$i;print $i}' > ${WRKDIR}/temp_all_samples.counts.txt
head -1 ${WRKDIR}/temp_all_samples.counts.txt > ${WRKDIR}/temp_header.txt
awk '{print $1}' ${WRKDIR}/final.out | grep -f - ${WRKDIR}/temp_all_samples.counts.txt > ${WRKDIR}/temp_retrocopy.txt
awk '{print $6}' ${WRKDIR}/final.out | grep -f - ${WRKDIR}/temp_all_samples.counts.txt > ${WRKDIR}/temp_parent.txt
cat ${WRKDIR}/temp_header.txt ${WRKDIR}/temp_retrocopy.txt ${WRKDIR}/temp_parent.txt > ${WRKDIR}/all_samples.counts.txt

rm -rf ${WRKDIR}/hisat2
rm -rf ${WRKDIR}/ballgown
rm ${WRKDIR}/temp*

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee -a $LOGFILE

##bash ~/Music/RetroFinder/bin/RNASeq.sh ./RANSeq/ GCF_genomic.fna retro.gtf 6 /home/wzy/biotools/hisat2-2.1.0/hisat2 /home/wzy/biotools/stringtie-1.3.3b.Linux_x86_64/stringtie /usr/bin/samtools /home/wzy/biotools/hisat2-2.1.0/hisat2-build ~/Music/Trichoplusia_Ni
