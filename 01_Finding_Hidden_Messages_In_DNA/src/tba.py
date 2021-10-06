
from typing import Dict, List


class TBA:

    """tba"""

    def __init__(self) -> None:
        pass

    def pattern_count(self, text: str, pattern: str) -> int:

        count = 0
        pattern_len = len(pattern)
        for i in range(0, len(text) - len(pattern) + 1):
            j = i + pattern_len
            if text[i:j] == pattern:
                count += 1

        return count

    def frequency_table(self, text: str, k: int) -> Dict[str, int]:
        """
        Function takes a string text and an integer k as input and returns their
        frequency table as a map of string keys to integer values.
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
        frequent_patterns = []
        frequent_map = self.frequency_table(text, k)

        max_frequency = max(frequent_map.values())
        for k, v in frequent_map.items():
            if frequent_map[k] == max_frequency:
                frequent_patterns.append(k)

        return frequent_patterns
