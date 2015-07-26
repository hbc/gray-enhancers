"""
the supplied enhancer sequences have two non-enhancer sequences on
the end, this cuts those off.

>-249753741_hs_Endogenous_chr1:6052584^WCGCGTYat-59^GATTGGYat-29^GATTGGYat14^YGCGGCKSat23
ACTGGCCGCTTCACTGaactccccagcagcctgtacgtttagtcctacccgggcCCGCCGCAggGATTGGCaccgcgagcgtttcgcgtcgggagctgaacccgagaGATTGGCaggcgccgggactgccgctgtcaGACGCGAccgcccaagaCACTGCGGCTCCTCA

to

>-249753741_hs_Endogenous_chr1:6052584^WCGCGTYat-59^GATTGGYat-29^GATTGGYat14^YGCGGCKSat23
aactccccagcagcctgtacgtttagtcctacccgggcCCGCCGCAggGATTGGCaccgcgagcgtttcgcgtcgggagctgaacccgagaGATTGGCaggcgccgggactgccgctgtcaGACGCGAccgcccaaga

this keeps only unique sequences
"""

from __future__ import print_function
import sys
import itertools

seen = set()

with open(sys.argv[1], "r") as in_handle:
    for name, seq in itertools.izip_longest(*[in_handle]*2):
        llen = len(seq)
        piece = seq[16:(llen-16)].strip()
        if piece not in seen:
            seen.update([piece])
            print(name.strip(), file=sys.stdout)
            print(piece, file=sys.stdout)
        else:
            print("%s sequence already seen, skipping %s." % (seq, name),
                  file=sys.stderr)
