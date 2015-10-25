## assign barcodes to enhancers
These are the scripts to look at enhancer + barcode statistics from the MPRA
experiment.

These scripts take the tiled enhancer sequences and creates a BWA index of the
tiles. Then they align the barcoded reads to the tiled sequences and determine
which barcode corresponds to which tile by determining which tile the barcode
aligns best to. Barcodes can have errors, so these scripts also do some
simple disabiguation of erroneous barcodes. They also drop barcodes that are
ambiguous.

The clean_read.py part of these scripts is highly dependent on the structure of
your data. What this does is remove barcode and some other static sequences from
the read, and places the barcode in the read name. This is so after alignment,
we can figure out what barcode went with the sequence. If you can provide FASTQ
files that have the barcode in the read name you can skip this step and run the
rest of the analysis.

Your reads should have the format `barcode-name` as the read name if you want
to skip the cleaning step.

## cleaning
Raw reads are expected to have the format:
```
-11mer barcode +
-12mer constant restriction site (TCTAGAGGTACC) +
-88bp enhancer +
-33bp constant seq (CAGTGAAGCGGCCAGTGATCGGAAGAGCACAC ) +
-6bp Truseq index (GTCTGA  or CGATGT) +
-41bp P7 sequence (ACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTT)
```

But if you can provide reads with names `barcode-name` you can skip this
step.


## How to run
1) Create a FASTA file of the enhancer sequences, enhancers.fa. The names of the
   sequences should not have any spaces in them.
2) index the sequence with bwa: bwa index enhancers.fa
3) stick barcode in the read name with: python ../code/clean_read.py sequences.fq
3) align sequences with bwa mem: bwa mem -t number-of-threads enhancers.fa sequences.enhancers.fq > alignments.sam
4) create a table of the barcode and best enhancer sequence alignment from the BAM file:
python bam_to_table.py alignments.sam
5) load into analysis_fixed.Rmd

## example analysis
```bash
infile=data/TN05_S1_L001_R1_001.fastq
cleaned=data/TN05_S1.cleaned.fastq
prefix=TN05
out_dir=new-data
enhancers=metadata/TN03_newstrategy_trimmed.fa

mkdir -p out-dir
python ../code/clean_read.py $infile > $cleaned
bwa mem -t 6 $enhancers $cleaned > $out_dir/$prefix.sam
python ../code/bam_to_table.py $out_dir/$prefix.sam
```

You can use this table to do a more in depth analysis.
The file analysis_fixed.Rmd has an example,

then run analysis_fixed.Rmd on the file in $out_dir/$prefix.tsv file.

## help
If you'd like to use this idea to do a similar experiment and these aren't
working for you post an issue and we'll work with you to get something
that will work with your data.
