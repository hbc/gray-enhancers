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

here is an example script

infile=data/TN05_S1_L001_R1_001.fastq
prefix=TN05
out_dir=new-data
enhancers=metadata/TN03_newstrategy_trimmed.fa

mkdir -p out-dir
bwa mem -t 6 $enhancers $infile > $out_dir/$prefix.sam
python ../code/single_bam_to_table.py $out_dir/$prefix.sam
then run analysis_fixed.Rmd on the file in $out_dir/$prefix.tsv file.
