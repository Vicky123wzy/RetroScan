# RetroScan

## Overview

RetroScan is an easy-to-use tool for retrocopy identification that integrates a series of bioinformatics tools (LAST, BEDtools, ClustalW2, KaKs_Calculator, HISAT2, StringTie, SAMtools and Shiny) and scripts. It scans retrocopies based on alignments between protein-coding genes and whole-genome sequences. This tool can also analyze heterosense substitution and synonymous substitution, compare gene structure between parental genes and retrocopies, and calculate corresponding expression values. Moreover, RetroScan has a user-friendly visualization interface that provides overall statistical information, a retrocopy structure diagram, the nonsynonymous/synonymous substitution (Ka/Ks) ratio distribution and the fragments per kilobase per million (FPKM) heatmap using the Shiny package in R.

It can be downloaded on https://github.com/Vicky123wzy/RetroScan.

![Aaron Swartz](https://github.com/Vicky123wzy/RetroScan/raw/main/pic/pipeline.png)

### Getting Started
#### Installation
##### conda:
```
conda create -n retroscan python=3.8 
conda activate retroscan
conda install -c biocanda -c biovicky retroscan
```
##### manully:
```
git clone https://github.com/Vicky123wzy/RetroScan.git
```
The retroscan.py is located in the /RetroScan/bin.


The RetroScan pipeline requires the following dependencies :
- `gffread == version 1.03.73` – https://github.com/gpertea/gffread
- `LAST == 1186`  -  http://last.cbrc.jp/
- `bedtools == v2.29.2`  -  http://last.cbrc.jp/
- `seqkit == 0.14.0`  -  https://bioinf.shenwei.me/seqkit/
- `diamond == v0.9.18.119`  -  http://github.com/bbuchfink/diamond
- `clustalw2 == 2.1`  -  http://www.clustal.org/clustal2/
- `KaKs_Calculator == 2.0`  -  https://sourceforge.net/projects/kakscalculator2
- `hisat2 == 2.2.1`  -  for RNA-Seq only (http://ccb.jhu.edu/software/hisat2)
- `stringtie == 2.1.4`  -  for RNA-Seq only (http://ccb.jhu.edu/software/stringtie)
- `samtools == 1.9`  -  for RNA-Seq only (http://samtools.sourceforge.net/)
- `Python >= 3.6`  -  require `panda` `matplotlib` `numpy`
- `Perl == v5.26.2`


### USAGE
example:
The example data is in the /RetroScan/example, or can be downloaded from https://github.com/Vicky123wzy/RetroScan/tree/main/example.

No RNA-Seq data:
```
retroscan.py test.fa.zip test.gff3.zip Output_dir
```
RNA-Seq data:
```
retroscan.py test.fa.zip test.gff3.zip Output_dir --RNASeq_path rnaseq_file
```


Other parameters:
```
retroscan.py [-h] [--identity [int]] [--coverage_rate [int]]
                    [--coverage_len [int]] [--intron_loss_num [int]]
                    [--gaplen [int]] [--deduplication_identity_cutoff [int]]
                    [--deduplication_coverage_cutoff [int]]
                    [--deduplication_similar_cutoff [int]]
                    [--kaksmethod [{NG,LWL,LPB,MLWL,MLPB,GY,YN,MYN,MS,MA,GNG,GLWL,GLPB,GMLWL,GMLPB,GYN,GMYN}]]
                    [--thread [int]] [--gffread [GFFREAD]] [--lastdb [LASTDB]]
                    [--lastal [LASTAL]] [--bedtools [BEDTOOLS]]
                    [--seqkit [SEQKIT]] [--diamond [DIAMOND]]
                    [--clustalw2 [CLUSTALW2]] [--RNASeq_path [RNASEQ_PATH]]
                    [--hisat2 [HISAT2]] [--hisat2build [HISAT2BUILD]]
                    [--stringtie [STRINGTIE]] [--samtools [SAMTOOLS]]
                    [--version]
                    genome_file gff_file Output_dir

positional arguments:
  genome_file                           - Genome sequences fasta file, ending with fasta, fa,fas, faa, fna, tar.gz or zip
  gff_file                              - Gff file, ending with gff, gff3, tar.gz or zip
  Output_dir                            - Path to the directory where all intermediate and final results will be stored

optional arguments:
  -h, --help                            - show this help message and exit
  --identity [int]                      - Expectation identity threshold(percentage) for saving hits [50]
  --coverage_rate [int]                 - Expectation coverage rate(percentage) threshold for saving hits [50]
  --coverage_len [int]                  - Expectation length(aa) of coverage for saving hits [50]
  --intron_loss_num [int]               - Expectation the number of intron-lost for saving hits [2]
  --gaplen [int]                        - Expectation length(bp) of gap for merging hits [40]
  --deduplication_identity_cutoff [int] - The identity cutoff (percentage) of retrocopies duplication sequences aligned to the genome [80]
  --deduplication_coverage_cutoff [int] - The coverage cutoff (percentage) of retrocopies duplication sequences aligned to the genome [80]
  --deduplication_similar_cutoff [int]  - The cutoff number of retrocopies duplication sequences aligned to the genome [10]
  --kaksmethod                          - Methods for estimating Ka and Ks and theirs references, choices=['NG', 'LWL', 'LPB', 'MLWL', 'MLPB', 'GY', 'YN', 'MYN', 'MS', 'MA', 'GNG', 'GLWL', 'GLPB', 'GMLWL', 'GMLPB', 'GYN', 'GMYN']
  --thread [int]                        - Number of threads (CPUs) to use [1]
  --gffread [GFFREAD]                   - Path to the gffread, https://github.com/gpertea/gffread (Default is "gffread").
  --lastdb [LASTDB]                     - Path to the lastdb, http://last.cbrc.jp/ (Default is "lastdb").
  --lastal [LASTAL]                     - Path to the lastal, http://last.cbrc.jp/ (Default is "lastal").
  --bedtools [BEDTOOLS]                 - Path to the bedtools, https://bedtools.readthedocs.io/en/latest/index.html (Default is "bedtools").
  --seqkit [SEQKIT]                     - Path to the seqkit, https://bioinf.shenwei.me/seqkit/ (Default is "seqkit").
  --diamond [DIAMOND]                   - Path to the diamond, http://github.com/bbuchfink/diamond (Default is "diamond").
  --clustalw2 [CLUSTALW2]               - Path to the clustalw2, http://www.clustal.org/clustal2/ (Default is "clustalw2").
  --RNASeq_path [RNASEQ_PATH]           - Path to the RNASeq file.
  --hisat2 [HISAT2]                     - Path to the hisat2,  https://daehwankimlab.github.io/hisat2/ (Default is "hisat2").
  --hisat2build [HISAT2BUILD]           - Path to the hisat2-build, https://daehwankimlab.github.io/hisat2/ (Default is "hisat2-build").
  --stringtie [STRINGTIE]               - Path to the stringtie, http://ccb.jhu.edu/software/stringtie/ (Default is "stringtie ").
  --samtools [SAMTOOLS]                 - Path to the samtools, http://www.htslib.org/ (Default is "samtools").
  --version                             - show program's version number and exit
```

### Input
If only two files: genome sequences file (fasta format) and corresponding annotation file (gff format)  are input into RetroScan, it can identify the detailed information of retrocopies and parental genes on the genome. 

The user needs to provide the RNA-Seq data to obtain the expression information of retrocopy. The folder structure is shown below:
```
 RANSeq file
    ├── sample1.fastq.gz
    ├── sample2_1.fastq.gz
    └── sample2_2.fastq.gz
```

### Output
The pipeline generates a number of files. The types of files are listed below. 

a) 	final.out
- 	the detailed information of retrocopies and parental genes, including: Retrocopy_ID, Retro_chr, Retro_strand, Retro_start, Retro_end, Parental_gene_ID, Parent_chr, Parent_start, Parent_end, Pro_start, Pro_end, Identity, Coverage, Lost_intron, Description, Ka, Ks, Ka/Ks, Host_gene_ID, Retro_sequence, Retro_protein. They are seperated by "\t".

b) 	gene_trans_info.txt
-	the detailed information of all protein-coding gene, including: gene, mRNA, mRNA_chr, mRNA_strand, mRNA_start, mRNA_end, protein_length, cds_num, cds_site, exon_site, pep_site. They are seperated by "\t".

c) 	protein.fasta
- 	protein sequences of the longest transcript of each gene. 

d) 	cds.fasta
- 	all cds sequences of all transcripts.

e) 	all_samples.counts.txt
- 	FPKM values for retrocopies and parental genes of each samples.

f) 	RetroScan.log
- 	the log file for all pipeline.


#### visualization
The online visualization website is: https://bioinfovicky.shinyapps.io/retroscan-app/.

User can upload the results files and download the figures.

Users can also download the visualization webpage (https://github.com/Vicky123wzy/RetroScan/tree/main/RetroScan-app) to compile and use in local R.
```
library(shiny)
library(dplyr)
library(stringr)
library(png)
library(shinyjs)
library(DT)
library(visNetwork)
library(rintrojs)
library(DT)
library(shinydashboard)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(UpSetR)
library(Biostrings)
library(muscle)
library(ape)
library(ggmsa)
library(patchwork)
library(pheatmap)
library(colourpicker)

shiny::runApp('RetroScan-app')
```


### Citations and licensing
If you use this code or the resulting assemblies, please cite the following  paper:

RetroScan: an easy-to-use pipeline for retrocopy annotation and visualization. 


#### This software is distributed under The MIT License (MIT).
