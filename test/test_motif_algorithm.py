import unittest
import numpy as np
from src.motif.motif import MotifAlgorithm


class TestMotifAlgorithm(unittest.TestCase):

    def __init__(self, methodName: str = ...) -> None:
        super().__init__(methodName=methodName)

        self.uut = MotifAlgorithm()
        self.dataPath = "data/replication"

    def read_data(self, file_path: str, read_type: str, file_count: int) -> list:
        data = []
        for i in range(1, file_count+1):
            with open(file_path+str(i)+".txt") as f:
                if read_type == "int":
                    line = f.read().splitlines()
                    data.append(list(map(int, line)))
                else:
                    data.append(f.read().splitlines())

        return data

    def test_motif_enumeration(self):
        inputs = self.read_data(self.dataPath+"/motif_enum_inputs/input_", "str", 7)
        expected = self.read_data(self.dataPath+"/motif_enum_outputs/output_", "str", 7)

        actual = []
        for input in inputs:
            kd = input[0].split(' ')
            genomes = list(input[1].split(' '))
            actual.append(self.uut.motif_enumeration(genomes, int(kd[0]), int(kd[1])))

        for e1, a1 in zip(expected, actual):
            if not e1:
                self.assertEqual(e1, a1)
            else:
                self.assertCountEqual(e1[0].split(' '), a1)

    def test_distance_between_pattern_and_strings(self):
        inputs = self.read_data(self.dataPath+"/distance_between_pattern_and_strings_inputs/input_", "str", 6)
        expected = self.read_data(self.dataPath+"/distance_between_pattern_and_strings_outputs/output_", "int", 6)

        actual = []
        for input in inputs:
            dna = list(input[1].split(' '))
            actual.append([self.uut.distance_between_pattern_and_strings(input[0], dna)])

        self.assertEqual(expected, actual)

    def test_median_string(self):
        inputs = self.read_data(self.dataPath+"/median_string_inputs/input_", "str", 5)
        expected = self.read_data(self.dataPath+"/median_string_outputs/output_", "str", 5)

        actual = []
        for input in inputs:
            dna = list(input[1].split(' '))
            actual.append(self.uut.median_string(dna, int(input[0])))

        for e1, a1 in zip(expected, actual):
            self.assertIn(e1[0], a1)

    def test_greedy_motif_search(self):
        inputs = self.read_data(self.dataPath+"/greedy_motif_search_inputs/input_", "str",6)
        expected = self.read_data(self.dataPath+"/greedy_motif_search_outputs/output_", "str", 6)

        actual = []
        for input in inputs:
            k = int(input[0].split(' ')[0])
            dna = list(input[1].split(' '))
            actual.append(self.uut.greedy_motif_search(dna, k))

        for e1, a1 in zip(expected, actual):
            self.assertEqual(list(e1[0].split(' ')), a1)

    def test_greedy_motif_pseudocount_search(self):
        pass

if __name__ == '__main__':
    unittest.main()
