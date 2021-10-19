
from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))
print(sys.path)


import math
import numpy as np
from collections import Counter
from src.replication.genome_replication_algorithm import GenomeReplicationAlgorithm
from typing import List



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
        """Score

        The total number of unpopular nucleotide in the motif column matrix
        """
        total = motifs.shape[0]

        score = 0
        for col in motifs.T:
            score += total - Counter(col).most_common(1)[0][1]

        return score

    def count(self, motifs: np.ndarray) -> np.ndarray:
        """Count

        The total number of nucleotide in the motif column matrix
        """
        res = np.zeros(shape=(motifs.shape[1], 4), dtype="int")
        for idx, col in enumerate(motifs.T):
            counts = Counter({'A': 0, 'C': 0, 'G': 0, 'T': 0})
            counts.update(col)
            res[idx] = [counts['A'], counts['C'], counts['G'], counts['T']]

        return res.T

    def profile(self, motifs: np.ndarray) -> np.ndarray:
        """Profile

        The probability of each nucleotide in the motif column matrix
        """
        dividend = motifs.shape[0]
        counts = self.count(motifs)
        return counts/dividend

    # TODO: return all combinations of tied breaker
    def consensus(self, motifs: np.ndarray) -> str:
        """Consensus String

        The most popular letters in each column of the motif matrix
        """
        pprofile = self.profile(motifs).T

        consensus_motif = []
        for p in pprofile:
            idx = np.argmax(p)
            if idx == 0:
                consensus_motif.append('A')
            elif idx == 1:
                consensus_motif.append('C')
            elif idx == 2:
                consensus_motif.append('G')
            elif idx == 3:
                consensus_motif.append('T')

        return ''.join(consensus_motif)

    def entropy(self, motifs: np.ndarray) -> float:
        """Entropy

        Total uncertainty of a probability distribution in the motif matrix
        """
        pprofile = self.profile(motifs).T

        total_entropy = 0
        for probs in pprofile:
            sub_entropy = 0
            for val in probs:
                if val != 0:
                    sub_entropy += val * math.log2(val)
            total_entropy -= sub_entropy

        return total_entropy

    def distance_between_pattern_and_strings(self, pattern: str, dna: List[str]) -> int:
        """d(pattern, dna)

        The sum of distances between Pattern and each string in Dna = {Dna1, ..., Dnat}
        """
        k = len(pattern)
        min_distance = 0

        for s in dna:
            min_hd = sys.maxsize
            for i in range(0, len(s) - k + 1):
                hd = self.gra.hamming_distance(pattern, s[i:i+k])
                if min_hd > hd:
                    min_hd = hd
            min_distance += min_hd

        return min_distance

    def median_string(self, dna: str, k: int) -> List[str]:
        """Median String


        Input:  An integer k, followed by a space-separated collection of
                strings Dna.
        Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all possible
                choices of k-mers.
        """
        distance = sys.maxsize
        patterns = self.gra.neighbors("A"*k, k)
        median = []
        for pattern in patterns:
            d = self.distance_between_pattern_and_strings(pattern, dna)
            if distance > d:
                distance = d
                median.append(pattern)

        return median


