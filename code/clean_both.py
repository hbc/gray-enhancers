import os
from argparse import ArgumentParser
from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists
from Bio import SeqIO
from itertools import izip

DEFAULT_FORMAT = "fastq-sanger"

def partition_read1(fastq_file_1, fastq_file_2):
    base_1, ext_1 = os.path.splitext(fastq_file_1)
    base_2, ext_2 = os.path.splitext(fastq_file_2)
    read_number = 0
    in_handle_1 = SeqIO.parse(fastq_file_1, "fastq-sanger")
    in_handle_2 = SeqIO.parse(fastq_file_2, "fastq-sanger")
    enhancer_fq_1 = base_1 + ".enhancer" + ext_1
    enhancer_fq_2 = base_2 + ".enhancer" + ext_2
    out_files = [enhancer_fq_1, enhancer_fq_2]
    if all(map(file_exists, out_files)):
        return out_files

    skipped_too_short = 0
    skipped_barcode_not_match = 0
    skipped_restriction_not_TCT = 0
    skipped_not_TGA = 0
    read_used = 0

    with file_transaction(out_files) as tx_out_files:
        tx_fq_out_1 = tx_out_files[0]
        tx_fq_out_2 = tx_out_files[1]
        with open(tx_fq_out_1, "w") as fq_handle_1, open(tx_fq_out_2, "w") as fq_handle_2:
            for read1, read2 in izip(in_handle_1, in_handle_2):
                read_number += 1
                if len(read1) < 187 or len(read2) < 187:
                    skipped_too_short += 1
                    continue
                barcode_1 = read1.seq[0:16]
                restriction = read1.seq[16:28]
                if not restriction.startswith("TCT"):
                    skipped_restriction_not_TCT += 1
                    continue
                adapter = read1.seq[28:44]
                if not adapter.startswith("TGA"):
                    skipped_not_TGA += 1
                    continue
		read2_desc = str(read2.description)
                read2 = read2.reverse_complement()
                barcode_2 = read2.seq[0:16]
                if str(barcode_1) != str(barcode_2):
                    skipped_barcode_not_match += 1
                    continue
                read_used += 1
                read1.id = "{barcode_1}-{read_number}".format(**locals())
                fq_handle_1.write(str(read1[44:187].__format__("fastq-sanger")))
                read2.id = "{barcode_2}-{read_number}".format(**locals())
		read2.description = read2_desc
                fq_handle_2.write(str(read2[44:187].__format__("fastq-sanger")))
    print "Processed %d total reads." % read_number
    print "Skipped %d reads because they were too short." % skipped_too_short
    print "Skipped %d reads because the R1 and R2 barcode didn't match.." % skipped_barcode_not_match
    print "Skipped %d reads because the restriction site didn't start with TCT." % skipped_restriction_not_TCT
    print "Skipped %d reads beacuse the adapter did not start with TGA." % skipped_not_TGA
    print "Kept %d reads." % read_used
    return enhancer_fq_1, enhancer_fq_2

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("fastq1", help="Read 1 of lane to parse.")
    parser.add_argument("fastq2", help="Read 2 of lane to parse.")

    args = parser.parse_args()

    partition_read1(args.fastq1, args.fastq2)

