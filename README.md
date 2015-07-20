# gray-enhancers

Scripts to look at enhancer + barcode statistics from the MPRA experiment

http://htmlpreview.github.io/?https://github.com/hbc/gray-enhancers/blob/master/code/analysis.html


## How to run
1) Create a FASTA file of the enhancer sequences, enhancers.fa
2) index the sequence with bwa: bwa index enhancers.fa
3) align sequences with bwa mem: bwa mem -t number-of-threads enhancers.fa sequences.fq > alignments.sam
4) create a table of the barcode and best enhancer sequence alignment from the BAM file:
python single_bam_to_table.py alignments.sam
5) load into analysis_fixed.Rmd

## assumptions of the code
this assumes the enhancer sequences are named in a certain way
single_bam_to_table.py does

here are all the steps you need to do once you have the indexed enhancers to align to

#!/bin/sh
#BSUB -u rory.kirchner@gmail.com
#BSUB -J bwa-aln
#BSUB -o bwa-aln.out
#BSUB -e bwa-aln.err
#BSUB -W 5:00
#BSUB -n 6
#BSUB -q mcore
#bwa mem -t 6 metadata/TN03_newstrategy_trimmed.fa data/TN05_S1_L001_R1_001.fastq > align/TN05_L001_R1_001.sam
#samtools view -bS align/TN05_L001_R1_001.sam > align/TN05_L001_R1_001.bam
#python ../code/single_bam_to_table.py align/TN05_L001_R1_001.sam > align/TN05_L001_R1.table
samtools sort align/TN05_L001_R1_001.bam TN05.sorted

then run analysis_fixed.Rmd on it
