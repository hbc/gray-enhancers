from argparse import ArgumentParser
from bcbio.bam import open_samfile
import pysam


OUT_HEADER = ["id", "barcode", "eid", "species", "subtype", "mismatch", "mapq", "as", "xs", 
              "clipped", "insertions", "deletions", "matched", "other"]

def parse_cigar_tuples(tuples):
    d = {"insertions": 0,
         "deletions": 0,
         "clipped": 0,
         "matched": 0,
         "other": 0}
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

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("BAM", help="BAMfile")
    args = parser.parse_args()

    print "\t".join(OUT_HEADER)

    with open_samfile(args.BAM) as in_file:
        for read in in_file:
            barcode, number = read.qname.split("-")
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
            out_line = map(str, [number, barcode, eid, species, subtype, mismatch, mapq, AS, XS, 
                                 clipped, insertions, deletions, matched, other])
            print "\t".join(out_line)
