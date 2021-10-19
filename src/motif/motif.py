from collections import Counter
from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))
print(sys.path)

from typing import List
from src.replication.genome_replication_algorithm import GenomeReplicationAlgorithm
import numpy as np

class MotifAlgorithm:

    def __init__(self) -> None:
        self.gra = GenomeReplicationAlgorithm()

    def motif_enumeration(self, dna_list: List[str], k: int, d: int) -> List[str]:
        """ Motif Enumeration

        Given a collection of strings DNA and an integer d, a k-mer is a
        (k,d)-motif if it appears in every string from Dna with at most d
        mismatches. Find all (k, d)-motifs in DNA.
        """

        patterns = []
        for i in range(0, len(dna_list[0]) - k + 1):
            patterns.append(dna_list[0][i:i+k])


        motifs = set()
        # for each k-mer pattern in the first string in Dna
        for pattern in patterns:
            # for each k-mer pattern’ differing from Pattern by at most d mismatches
            for pattern_d in self.gra.neighbors(pattern, d):
                # if pattern' appears in each string from DNA with at most d mismatches
                contains_d = [self.gra.approximate_pattern_count(pattern_d, dna, d) for dna in dna_list]
                if all([count > 0 for count in contains_d]):
                    motifs.add(pattern_d)

        return list(motifs)

    def score(self, motifs: np.ndarray) -> int:

        total = motifs.shape[0]
        motifs_t = motifs.T

        score = 0
        for col in motifs_t:
            score += total - Counter(col).most_common(1)[0][1]

        return score

    def count(self, motifs: np.ndarray) -> np.ndarray:

        motifs_t = motifs.T

        res = np.zeros(shape=(motifs.shape[1], 4), dtype="int")


        for idx, col in enumerate(motifs_t):
            counts = Counter({'A':0, 'C':0, 'G':0, 'T':0})
            counts.update(col)
            res[idx] = [counts['A'], counts['C'], counts['G'], counts['T']]


        return res.T


uut = MotifAlgorithm()
data = np.loadtxt('data/replication/motifs.txt', dtype="str")
print(data.shape[0])
print(uut.score(data))
print(uut.count(data))