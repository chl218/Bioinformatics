
from typing import Dict, List


class PatternFindingAlgorithm:

    """ Finding Hidden Messages in DNA
    """

    def __init__(self) -> None:
        pass

    def pattern_count(self, text: str, pattern: str) -> int:
        """ Pattern Count

        The number of times that a k-mer pattern appears as a substring of text.
        """
        count = 0
        pattern_len = len(pattern)
        for i in range(0, len(text) - len(pattern) + 1):
            j = i + pattern_len
            if text[i:j] == pattern:
                count += 1

        return count

    def frequent_words(self, text: str, k: int) -> List[str]:
        """ Frequent Words

        Finding the most frequent k-mers in a string text
        """

        frequent_patterns = set()
        counts = []
        n = len(text) - k + 1

        for i in range(0, n):
            counts.append(self.pattern_count(text, text[i:i+k]))

        max_count = max(counts)
        for i in range(0, n):
            if counts[i] == max_count:
                frequent_patterns.add(text[i:i+k])

        return list(frequent_patterns)

    def frequency_table(self, text: str, k: int) -> Dict[str, int]:
        """ Frequency Table

        Find the frequency of k-mers in text
        """

        freq_map = {}
        n = len(text) - k + 1

        for i in range(0, n):
            pattern = text[i:i+k]
            if pattern in freq_map:
                freq_map[pattern] += 1
            else:
                freq_map[pattern] = 1

        return freq_map

    def better_frequent_words(self, text: str, k: int) -> List[str]:
        """ Frequent Words (Optimized)

        Finding the most frequent k-mers in a string text
        """

        frequent_patterns = []
        frequent_map = self.frequency_table(text, k)

        max_frequency = max(frequent_map.values())
        for k, v in frequent_map.items():
            if frequent_map[k] == max_frequency:
                frequent_patterns.append(k)

        return frequent_patterns

    def reverse_complement(self, pattern: str) -> str:
        """Reverse Complement

        Find the reverse complement of a DNA string
        """
        rc = []
        for s in reversed(pattern):
            if s == "A":
                rc.append("T")
            elif s == "T":
                rc.append("A")
            elif s == "G":
                rc.append("C")
            elif s == "C":
                rc.append("G")

        return ''.join(rc)

    def pattern_match(self, pattern: str, genome: str) -> List[int]:
        """Pattern Matching

        Returns a list of integers specifying all starting positions where
        pattern appears as a substring of genome.
        """
        idx = []

        n = len(genome) - len(pattern) + 1
        k = len(pattern)
        for i in range(0, n):
            if genome[i:i+k] == pattern:
                idx.append(i)

        return idx

    def find_clumps(self, text: str, k: int, L: int, t: int) -> List[str]:
        """Find Clumps

        Given integers L and t, a k-mer Pattern forms an (L, t)-clump inside a
        string Genome if there is an interval of Genome of length L in which
        this k-mer appears at least t times.
        """

        patterns = set()
        n = len(text) - L + 1
        nn = len(text)
        for i in range(0, n):
            print("%.3f" % (i/nn))
            window = text[i:i+L]
            freq_map = self.frequency_table(window, k)

            for key, val in freq_map.items():
                if val >= t:
                    patterns.add(key)

        return list(patterns)
