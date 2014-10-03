#!/bin/bash
#bwa index metadata/Array3_12k_2014_05_12_merged.fa
mkdir align
bwa bwasw metadata/Array3_12k_2014_05_12_merged.fa data/TN03_S1_L001_R1_001.subset.enhancer.fastq | samtools view -Sbh - > align/TN03_S1_L001_R1_001.subset.enhancer.bam