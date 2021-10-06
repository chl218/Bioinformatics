
from typing import Dict, List, Set


class DNAPatternFinding:

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

        Takes a string text and an integer k as input and returns their frequency
        table as a map of string keys to integer values.
        """

        freq_map = {}
        n = len(text)

        for i in range(0, n-k+1):
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
