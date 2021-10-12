from typing import List
from src.replication.pattern_finding_algorithm import PatternFindingAlgorithm


class GenomeReplicationAlgorithm:

    def __init__(self) -> None:
        self.patternFinding = PatternFindingAlgorithm()

    def skew_GC(self, genome: str) -> List[int]:
        """Skew Diagram

        Skew diagram between G C nucleotides
        """

        diff = 0
        skew = [0]

        n = len(genome)
        for i in range(0, n):
            if genome[i] == 'G':
                diff += 1
                skew.append(diff)
            elif genome[i] == 'C':
                diff -= 1
                skew.append(diff)
            else:
                skew.append(diff)

        return skew

    def skew_min(self, genome: str) -> List[int]:
        """Minimum Skew

        Find the position(s) in a genome where the skew diagram attains a minimum.
        """

        skew = self.skew_GC(genome)
        mmin = min(skew)
        return [idx for idx, val in enumerate(skew) if val == mmin]

    def hamming_distance(self, p: str, q: str) -> int:
        """ Hamming Distance

        Find the number of mismatches between strings p and q
        """
        diff = 0
        for pi, qi in zip(p, q):
            if pi != qi:
                diff += 1

        return diff

    def approximate_pattern_match(self, pattern: str, text: str, d: int) -> List[int]:
        """ Approximate Pattern Match

        Find all approximate occurrences of a pattern in a string where pattern
        appears as a substring of text with at most d mismatches.
        """
        n = len(text) - len(pattern) + 1
        k = len(pattern)

        positions = []
        for i in range(0, n):
            p2 = text[i:i+k]
            if self.hamming_distance(pattern, p2) <= d:
                positions.append(i)

        return positions

    def approximate_pattern_count(self, pattern: str, text: str, d: int) -> int:
        """ Approximate Pattern Count

        Given strings Text and Pattern as well as an integer d, count the total
        number of occurrences of pattern in Ttext with at most d mismatches
        """

        n = len(text) - len(pattern) + 1
        k = len(pattern)

        counts = 0
        for i in range(0, n):
            p2 = text[i:i+k]
            if self.hamming_distance(pattern, p2) <= d:
                counts += 1

        return counts

    def neighbors(self, pattern: str, d: int) -> List[str]:
        """ Neighbors

        Generate d-neighborhood of all k-mers whose Hamming distance from
        pattern does not exceed d
        """
        if d == 0:
            return [pattern]
        if len(pattern) == 1:
            return ['A', 'C', 'G', 'T']

        neighbor = set()
        suffix_pattern = pattern[1:]
        suffix_neighbors = self.neighbors(suffix_pattern, d)
        for sn in suffix_neighbors:
            if self.hamming_distance(suffix_pattern, sn) < d:
                neighbor.add('A' + sn)
                neighbor.add('C' + sn)
                neighbor.add('G' + sn)
                neighbor.add('T' + sn)
            else:
                neighbor.add(pattern[0] + sn)

        return list(neighbor)

    def frequent_words_with_mismatches(self, text: str, k: int, d: int) -> List[str]:
        """ Most Frequent Words with Mismatch

        A most frequent k-mer with up to d mismatches
        """
        freq_map = {}

        n = len(text) - k + 1
        for i in range(0, n):
            pattern = text[i:i+k]
            neighborhood = self.neighbors(pattern, d)
            for neighbor in neighborhood:
                if neighbor in freq_map:
                    freq_map[neighbor] += 1
                else:
                    freq_map[neighbor] = 1

        most_freq = []
        max_freq = max(freq_map.values())
        for k, v in freq_map.items():
            if v == max_freq:
                most_freq.append(k)

        return most_freq

    def frequent_words_with_mismatches_and_reverse_complement(self, text: str, k: int, d: int) -> List[str]:
        """ Most Frequent Words with Mismatch and Reverse Complement

        Find the most frequent k-mers (with mismatches and reverse complements) in a string.
        """
        freq_map = {}

        n = len(text) - k + 1
        for i in range(0, n):
            pattern = text[i:i+k]
            pattern_rc = self.patternFinding.reverse_complement(pattern)

            neighborhood = self.neighbors(pattern, d) + self.neighbors(pattern_rc, d)
            for neighbor in neighborhood:
                if neighbor in freq_map:
                    freq_map[neighbor] += 1
                else:
                    freq_map[neighbor] = 1

        most_freq = []
        max_freq = max(freq_map.values())
        for k, v in freq_map.items():
            if v == max_freq:
                most_freq.append(k)

        return most_freq
