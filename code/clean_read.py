import os
from argparse import ArgumentParser
from Bio import SeqIO, pairwise2

def file_exists(fname):
    """Check if a file exists and is non-empty.
    """
    try:
        return fname and os.path.exists(fname) and os.path.getsize(fname) > 0
    except OSError:
        return False

def find_sequence_index(seq, adapter):
    """
    if the read is missing the adapter sequence, then drop it
    """
    # first alignment is best alignment
    alignment = pairwise2.align.localxs(adapter, seq, -1, -0.1)[0]
    adapter_idx = start_from_alignment(alignment, len(adapter))
    if not adapter_idx:
        return None
    return adapter_idx

def start_from_alignment(alignment, adapter_length):
    if not alignment:
        return None
    # if the alignment is poor, skip it (there is a gap or two mismatches)
    # 14 = two mismatches or 1 gap with our parameters
    if alignment[2] < (adapter_length - 2):
        return None
    return alignment[3]

def partition_read1(seq, args):
    left_adapter_idx = find_sequence_index(seq, args.left_adapter)
    right_adapter_idx = find_sequence_index(seq, args.right_adapter)
    restriction_idx = find_sequence_index(seq, args.restriction)
    if not all([left_adapter_idx, right_adapter_idx, restriction_idx]):
        return None, None
    if restriction_idx < 16:
        return None, None
    barcode = seq[restriction_idx-16:restriction_idx]
    return barcode, (left_adapter_idx, right_adapter_idx + len(args.right_adapter))

def partition_reads(args):
    fastq_file_1 = args.fastq1
    base_1, ext_1 = os.path.splitext(fastq_file_1)
    read_number = 0
    in_handle_1 = SeqIO.parse(fastq_file_1, "fastq-sanger")
    enhancer_fq_1 = base_1 + ".enhancer" + ext_1
    if file_exists(enhancer_fq_1):
        return enhancer_fq_1

    skipped_too_short = 0
    skipped_missing_constant_seqs = 0
    kept = 0

    with open(enhancer_fq_1, "w") as fq_handle_1:
        for read in in_handle_1:
            read_number += 1
            if read_number % 100 == 0:
                print "Processed %d reads." % read_number
            if len(read) < 187:
                skipped_too_short += 1
                continue
            barcode, idx = partition_read1(read.seq, args)
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
    parser.add_argument("--left-adapter", default="TGAGGAGCCGCAGTG")
    parser.add_argument("--right-adapter", default="CAGTGAAGCGGCCAG")
    parser.add_argument("--restriction", default="TCTAGAGGTACC")
    args = parser.parse_args()

    partition_reads(args)
