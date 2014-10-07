from argparse import ArgumentParser
from bcbio.bam import open_samfile
import pysam


OUT_HEADER = ["id", "barcode", "eid", "species", "subtype", "mismatch", "mapq"]


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
            out_line = map(str, [number, barcode, eid, species, subtype, mismatch, mapq])
            print "\t".join(out_line)
