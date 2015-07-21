from collections import Counter, defaultdict
from argparse import ArgumentParser
import os
import pysam
import toolz as tz
import sys


OUT_HEADER = ["id", "ename", "barcode", "eid", "species", "subtype", "mismatch",
              "mapq", "as", "xs", "clipped", "insertions", "deletions",
              "matched", "other", "seq"]

def open_samfile(in_file):
    if is_bam(in_file):
        return pysam.Samfile(in_file, "rb")
    elif is_sam(in_file):
        return pysam.Samfile(in_file, "r")
    else:
        raise IOError("in_file must be either a BAM file or SAM file. Is the "
                      "extension .sam or .bam?")

def is_bam(in_file):
    _, ext = os.path.splitext(in_file)
    if ext == ".bam":
        return True
    else:
        return False


def is_sam(in_file):
    _, ext = os.path.splitext(in_file)
    if ext == ".sam":
        return True
    else:
        return False

def parse_cigar_tuples(tuples):
    d = {"insertions": 0,
         "deletions": 0,
         "clipped": 0,
         "matched": 0,
         "other": 0}
    if not tuples:
        d["clipped"] = "NA"
        d["insertions"] = "NA"
        d["deletions"] = "NA"
        d["matched"] = "NA"
        d["other"] = "NA"
        return d
    for t in tuples:
        if t[0] == 4 or t[0] == 5:
            d["clipped"] += t[1]
        elif t[0] == 1:
            d["insertions"] += t[1]
        elif t[0] == 2:
            d["deletions"] += t[1]
        elif t[0] == 0:
            d["matched"] += t[1]
        else:
            d["other"] += t[1]
    return d

def is_clipped(read):
    cigar = parse_cigar_tuples(read.cigar)
    if cigar["clipped"] > 0:
        return True
    else:
        return False

def get_insertions(read):
    cigar = read.cigar
    insertions = []
    index = 0
    for t in cigar:
        if t[0] == 1:
            for x in range(t[1]):
                insertions.append(("I", index, read.seq[index]))
                index += 1
        elif t[0] == 2:
            for x in range(t[1]):
                index -= 1
        else:
            index += t[1]
    return insertions

def get_variants(read):
    variants = []
    md = partition_md(read)
    index = 0
    for x in md:
        if x.isdigit():
            index += int(x)
        elif x.startswith("^"):
            deleted = x[1:]
            for c in deleted:
                variants.append(("D", index, c))
                index += 1
        else:
            for c in x:
                variants.append(("M", index, c))
                index += 1
    insertions = get_insertions(read)
    variants += insertions
    return set(variants)

is_digit = lambda x: x.isdigit()

def partition_md(read):
    try:
        md = read.opt("MD")
    except KeyError:
        return None
    return ["".join(x) for x in tz.partitionby(is_digit, md)]

def diffs(read, mate):
    read_md = partition_md(read)
    mate_md = partition_md(mate)

def reconstruct_reference(read):
    index = 0
    ref = ""
    seq = read.seq
    partitions = partition_md(read)
    for p in partitions:
        if p.isdigit():
            ref = ref + seq[index:index + p]
            index = index + p

def get_errors(read, mate):
    read_var = get_variants(read)
    mate_var = get_variants(mate)
    all_var = read_var.union(mate_var)

    synthesis = read_var.intersection(mate_var)
    sequencing = all_var.difference(synthesis)
    return {"sequencing": sequencing, "synthesis": synthesis}

def bam_to_tab(in_file):
    out_file = os.path.splitext(in_file)[0] + ".tsv"
    # if os.path.exists(out_file):
    #     return out_file
    skipped_no_cigar = 0
    pairs_processed = 0
    skipped_no_mate = 0
    skipped_mismatch_position = 0
    skipped_clipped = 0
    skipped_unmapped = 0
    skipped_secondary = 0
    skipped_read_not_matched = 0
    synthesis_dict = defaultdict(Counter)
    sequencing_dict = defaultdict(Counter)
    seen = Counter()

    with open_samfile(in_file) as in_file, open(out_file, "w") as out_handle:
        print >>out_handle, "\t".join(OUT_HEADER)
        for read in in_file:

            try:
                ename = in_file.getrname(read.tid)
            except:
                ename = "unmapped"
            if read.is_secondary:
                skipped_secondary += 1
                continue

            barcode, number = read.qname.split("-")
            fragment_id = (barcode, ename)
            mapq = int(read.mapq)
            try:
                ename = in_file.getrname(read.tid)
                species = ename.split("_")[1]
                subtype = ename.split("_")[2]
                eid = ename.split("_")[0].replace("-", "")
            except ValueError:
                ename = "unmapped"
                species = "NA"
                subtype = "NA"
                eid = "NA"
            try:
                mismatch = read.opt("NM")
            except KeyError:
                mismatch = 0
            try:
                AS = read.opt("AS")
            except KeyError:
                AS = 0
            try:
                XS = read.opt("XS")
            except KeyError:
                XS = 0
            cigar = parse_cigar_tuples(read.cigar)
            clipped = cigar["clipped"]
            insertions = cigar["insertions"]
            deletions = cigar["deletions"]
            matched = cigar["matched"]
            other = cigar["other"]
            seq = str(read.seq)
            out_line = map(str, [number, ename, barcode, eid, species, subtype,
                                 mismatch, mapq, AS, XS, clipped, insertions,
                                 deletions, matched, other, seq])
            out_string = "\t".join(out_line)
            print >>out_handle, out_string
    print >>sys.stdout, "Total pairs processed: %d" % pairs_processed
    print >>sys.stdout, "Skipped due to having no mate: %d" % skipped_no_mate
    print >>sys.stdout, "Skipped due to mapping to different locations: %d" % skipped_mismatch_position
    print >>sys.stdout, "Skipped due to clipping of read: %d" % skipped_clipped
    print >>sys.stdout, "Skipped due to pair not matching: %d" % skipped_read_not_matched
    print >>sys.stdout, "Skipped the secondary alignment: %d" % skipped_secondary

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("BAM", help="BAMfile")
    args = parser.parse_args()

    bam_to_tab(args.BAM)
