
from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))


import math
import random
import numpy as np
import sys
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

    def score(self, motifs: List[str]) -> int:
        """Score

        The total number of unpopular nucleotide in the motif column matrix
        """
        total = len(motifs)
        column_vectors = [[sublist[idx] for sublist in motifs] for idx in range(len(motifs[0]))]

        score = 0
        for col in column_vectors:
            score += total - Counter(col).most_common(1)[0][1]

        return score

    def count(self, motifs: List[List[str]]) -> np.ndarray:
        """Count

        The total number of nucleotide in the motif column matrix
        """

        res = np.zeros(shape=(len(motifs[0]), 4), dtype="int")
        column_vectors = [[sublist[idx] for sublist in motifs] for idx in range(len(motifs[0]))]
        for idx, col in enumerate(column_vectors):
            counts = Counter({'A': 0, 'C': 0, 'G': 0, 'T': 0})
            counts.update(col)
            res[idx] = ([counts['A'], counts['C'], counts['G'], counts['T']])

        return res.T

    def profile(self, motifs: List[str]) -> dict:
        """Profile

        The probability of each nucleotide in the motif column matrix
        """

        dividend = len(motifs)
        probs = self.count(motifs)/dividend

        res = {'A': probs[0],
               'C': probs[1],
               'G': probs[2],
               'T': probs[3]}
        return res

    def pseudocount_profile(self, motifs: List[str]) -> dict:
        """ Profile with pseudocounts

        Laplace’s Rule of Succession in to improve unfair scoring
        """
        dividend = 2*len(motifs)
        pseudocount = 1 + self.count(motifs)
        probs = pseudocount/dividend

        res = {'A': probs[0],
               'C': probs[1],
               'G': probs[2],
               'T': probs[3]}
        return res

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
            if d < distance:
                distance = d
                median = [pattern]
            elif d == distance:
                median.append(pattern)

        return median

    def profile_most_probable_kmer(self, dna: str, k: int, profile_matrix: dict) -> str:
        """ Porilfe most probable k-mer

        Given a profile matrix Profile, we can evaluate the probability of every
        k-mer in a string Text and find a Profile-most probable k-mer in Text
        """
        p_max = 0
        kmer = dna[0:k]
        for i in range(0, len(dna) - k + 1):
            substr = dna[i:i+k]
            p = profile_matrix[substr[0]][0]
            for j in range(1, k):
                p *= profile_matrix[substr[j]][j]

            if p > p_max:
                p_max = p
                kmer = substr

        return kmer

    def greedy_motif_search(self, dna_list: List[str], k: int) -> List[str]:
        """ Greedy Motif Search

        Input:  Integers k and t, followed by a space-separated collection of
                strings Dna.
        Output: A collection of strings BestMotifs resulting from applying
                GreedyMotifSearch(Dna, k, t).
        """

        # BestMotifs ← motif matrix formed by first k-mers in each string from Dna
        best_motifs = []
        for dna in dna_list:
            best_motifs.append(dna[0:k])

        t = len(dna_list)
        n = len(dna_list[0]) - k + 1

        # for each k-mer Motif in the first string from Dna
        for j in range(0, n):
            motifs = [dna_list[0][j:j+k]]

            # form Profile from motifs Motif_1, ..., Motif_i-1
            for i in range(1, t):
                pprofile = self.profile(motifs)
                motifs.append(self.profile_most_probable_kmer(dna_list[i], k, pprofile))

            if self.score(motifs) < self.score(best_motifs):
                best_motifs = motifs

        return best_motifs

    def greedy_motif_pesudocount_search(self, dna_list: List[str], k: int) -> List[str]:
        """ Greedy motif search with pesudocount

        """
        best_motifs = []
        for dna in dna_list:
            best_motifs.append(dna[0:k])

        t = len(dna_list)
        n = len(dna_list[0]) - k + 1

        for j in range(0, n):
            motifs = [dna_list[0][j:j+k]]

            for i in range(1, t):
                #  Apply Laplace's Rule of Succession to form Profile from
                #  motifs Motif_1, ..., Motif_i-1
                pprofile = self.pseudocount_profile(motifs)
                motifs.append(self.profile_most_probable_kmer(dna_list[i], k, pprofile))

            if self.score(motifs) < self.score(best_motifs):
                best_motifs = motifs

        return best_motifs


    def motifs(self, profile: dict, dna_list: List[str]) -> List[str]:

        k = len(profile['A'])
        res = []
        for dna in dna_list:
            res.append(self.profile_most_probable_kmer(dna, k, profile))

        return res

    def randomized_motif_search(self, dna_list: List[str], k: int) -> List[str]:
        best_motifs = []
        for dna in dna_list:
            i = random.randint(0, len(dna)- k)
            best_motifs.append(dna[i:i+k])

        t = len(dna_list)
        while True:
            pprofile = self.pseudocount_profile(best_motifs)
            curr_motifs = self.motifs(pprofile, dna_list)

            if self.score(curr_motifs) < self.score(best_motifs):
                best_motifs = curr_motifs
            else:
                return best_motifs

    def randomized_motif_search_iteration(self, n: int, dna_list: List[str], k: int) -> List[str]:
        best_motif = self.randomized_motif_search(dna_list, k)
        for _ in range(n):
            curr_motif = self.randomized_motif_search(dna_list, k)
            if self.score(curr_motif) < self.score(best_motif):
                best_motif = curr_motif

        return best_motif

    def gibbs_sampler(self) -> List[str]:
        pass

    def rand_cdf(self, probs: List[float]) -> int:
        s = sum(probs)
        p_normalized = [p / s for p in probs]
        choices = list(range(1, len(probs)+1))
        return np.random.choice(choices, p=p_normalized)



uut = MotifAlgorithm()

p = [0.75, 0.50, 0.25, 0.1]
res = []
for i in range(100):
    res.append(uut.rand_cdf(p))

print(Counter(res))