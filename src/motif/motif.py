from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))
print(sys.path)

from typing import List
from src.replication.genome_replication_algorithm import GenomeReplicationAlgorithm

class MotifAlgorithm:

    def __init__(self) -> None:

        self.gra = GenomeReplicationAlgorithm()

    def motif_enumeration(self, dna_list: List[str], k: int, d: int) -> List[str]:
        """ Motif Enumeration

        Given a collection of strings DNA and an integer d, a k-mer is a
        (k,d)-motif if it appears in every string from Dna with at most d
        mismatches. Find all (k, d)-motifs in DNA.
        """

        motifs = set()
        dna = ''.join(dna_list)
        n = len(dna) - k + 1

        # for each k-mer Pattern in the first string in Dna
        for i in range(0, n):
            pattern = dna[i:i+k]
            neighbors = self.gra.neighbors(pattern, d)

            # for each k-mer Pattern’ differing from Pattern by at most d mismatches
            for pattern_d in neighbors:
                # if Pattern' appears in each string from Dna with at most d mismatches
                containsd = [self.gra.approximate_pattern_count(pattern_d, dna, d) for dna in dna_list]
                if all([count > 0 for count in containsd]):
                    motifs.add(pattern_d)

        return list(motifs)


