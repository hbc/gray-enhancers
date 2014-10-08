from argparse import ArgumentParser


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--trim", help="bases to trim from front and back of read",
                        default=16)
    parser.add_argument("fasta", help="FASTA file of enhancer sequences")
    args = parser.parse_args()

    header = ["eid", "species", "subtype", "sequence"]

    with open(args.fasta) as in_handle:
        for line in in_handle:
            if line.startswith(">"):
                ename = line.replace(">", "")
                species = ename.split("_")[1]
                subtype = ename.split("_")[2]
                eid = ename.split("_")[0].replace("-", "")
                continue

            seq = line[args.trim:-args.trim].lower()
            print "\t".join([eid, species, subtype, seq])
