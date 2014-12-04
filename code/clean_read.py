import os
import re
from argparse import ArgumentParser
from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists
from Bio import SeqIO, pairwise2
from itertools import izip, takewhile

DEFAULT_FORMAT = "fastq-sanger"
LEFT_ADAPTER_SEQUENCE = "TGAGGAGCCGCAGTG"
RIGHT_ADAPTER_SEQUENCE = "CAGTGAAGCGGCCAG"
RESTRICTION_SITE = "TCTAGAGGTACC"

def find_sequence_index(seq, adapter):
    """
    if the read is missing the adapter sequence, then drop it
    """
    # first alignment is best alignment
    alignment = pairwise2.align.localxs(adapter, seq, -1, -0.1)[0]
    adapter_idx = start_from_alignment(alignment)
    if not adapter_idx:
        return None
    return adapter_idx

def start_from_alignment(alignment):
    if not alignment:
        return None
    # if the alignment is poor, skip it (there is a gap or two mismatches)
    # 14 = two mismatches or 1 gap with our parameters
    if alignment[2] < 14:
        return None
    return alignment[3]

def partition_read1(seq):
    left_adapter_idx = find_sequence_index(seq, LEFT_ADAPTER_SEQUENCE)
    right_adapter_idx = find_sequence_index(seq, RIGHT_ADAPTER_SEQUENCE)
    restriction_idx = find_sequence_index(seq, RESTRICTION_SITE)
    if not all([left_adapter_idx, right_adapter_idx, restriction_idx]):
        return None, None
    if restriction_idx < 16:
        return None, None
    barcode = seq[restriction_idx-16:restriction_idx]
    return barcode, (left_adapter_idx, right_adapter_idx + len(RIGHT_ADAPTER_SEQUENCE))

def partition_reads(fastq_file_1):
    base_1, ext_1 = os.path.splitext(fastq_file_1)
    read_number = 0
    in_handle_1 = SeqIO.parse(fastq_file_1, "fastq-sanger")
    enhancer_fq_1 = base_1 + ".enhancer" + ext_1
    if file_exists(enhancer_fq_1):
        return enhancer_fq_1

    skipped_too_short = 0
    skipped_missing_constant_seqs = 0
    kept = 0

    with file_transaction(enhancer_fq_1) as tx_out_file:
        with open(tx_out_file, "w") as fq_handle_1:
            for read in in_handle_1:
                read_number += 1
                if len(read) < 187:
                    skipped_too_short += 1
                    continue
                barcode, idx = partition_read1(read.seq)
                if not barcode or not idx:
                    skipped_missing_constant_seqs += 1
                    continue
                kept += 1
                read.id = "{barcode}-{read_number}".format(**locals())
                fq_handle_1.write(str(read[idx[0]:idx[1]].__format__("fastq-sanger")))

    print "Reads proccessed: %s." % read_number
    print "Reads skipped for being too short: %s." % skipped_too_short
    print ("Reads skipped for not matching the anchor sequences: "
           "%s." % skipped_missing_constant_seqs)
    print "Reads kept: %s." % kept
    return enhancer_fq_1

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("fastq1", help="Read 1 of lane to parse.")

    args = parser.parse_args()

    partition_reads(args.fastq1)
