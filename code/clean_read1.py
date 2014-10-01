import os
from argparse import ArgumentParser
from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists
from Bio import SeqIO

DEFAULT_FORMAT = "fastq-sanger"

def partition_read1(fastq_file):
    base, ext = os.path.splitext(fastq_file)
    read_number = 0
    in_handle = SeqIO.parse(fastq_file, "fastq-sanger")
    enhancer_fq = base + ".enhancer" + ext
    enhancer_tsv = base + ".enhancer.tsv"
    out_files = [enhancer_fq, enhancer_tsv]
    if file_exists(enhancer_fq) and file_exists(enhancer_tsv):
        return enhancer_fq, enhancer_tsv

    with file_transaction(out_files) as tx_out_files:
        tx_fq_out = tx_out_files[0]
        tx_tsv_out = tx_out_files[1]
        with open(tx_fq_out, "w") as fq_handle, open(tx_tsv_out, "w") as tsv_handle:
            for read in in_handle:
                read_number += 1
                if len(read) < 187:
                    continue
                barcode = read.seq[0:16]
                restriction = read.seq[16:28]
                if not restriction.startswith("TCT"):
                    continue
                adapter = read.seq[28:44]
                if not adapter.startswith("TGA"):
                    continue
                enhancer = read.seq[44:187]
                out_line = "\t".join(map(str, [barcode, restriction, adapter, enhancer]))
                tsv_handle.write(out_line + "\n")
                read.id = "{barcode}-{read_number}".format(**locals())
                fq_handle.write(str(read[44:187].__format__("fastq-sanger")))
    return enhancer_fq, enhancer_tsv

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("fastq1", help="Read 1 of lane to parse.")

    args = parser.parse_args()

    partition_read1(args.fastq1)
