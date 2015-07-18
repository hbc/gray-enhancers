"""
the supplied enhancer sequences have two non-enhancer sequences on
the end, this cuts those off.

>-249753741_hs_Endogenous_chr1:6052584^WCGCGTYat-59^GATTGGYat-29^GATTGGYat14^YGCGGCKSat23
ACTGGCCGCTTCACTGaactccccagcagcctgtacgtttagtcctacccgggcCCGCCGCAggGATTGGCaccgcgagcgtttcgcgtcgggagctgaacccgagaGATTGGCaggcgccgggactgccgctgtcaGACGCGAccgcccaagaCACTGCGGCTCCTCA

to

>-249753741_hs_Endogenous_chr1:6052584^WCGCGTYat-59^GATTGGYat-29^GATTGGYat14^YGCGGCKSat23
aactccccagcagcctgtacgtttagtcctacccgggcCCGCCGCAggGATTGGCaccgcgagcgtttcgcgtcgggagctgaacccgagaGATTGGCaggcgccgggactgccgctgtcaGACGCGAccgcccaaga
"""

import sys

print sys.argv[1]
with open(sys.argv[1], "r") as in_handle:
    for line in in_handle:
        if line.startswith(">"):
            print line.strip()
        else:
            llen = len(line)
            print line[16:(llen-16)].strip()
