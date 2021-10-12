import unittest
import numpy as np
from src.replication.genome_replication_algorithm import GenomeReplicationAlgorithm


class TestGenomeReplicationAlgorithm(unittest.TestCase):

    def __init__(self, methodName: str = ...) -> None:
        super().__init__(methodName=methodName)

        self.uut = GenomeReplicationAlgorithm()
        self.dataPath = "data/replication"

    def read_data(self, file_path: str, read_type: str, file_count: int) -> list:
        data = []
        for i in range(1, file_count+1):
            res = np.loadtxt(file_path+str(i)+".txt", dtype=read_type)
            if res.size == 1:
                data.append([res.tolist()])
            else:
                data.append(res.tolist())
        return data

    def test_skew_min(self):
        inputs = self.read_data(self.dataPath+"/minimum_skew_inputs/input_", "str", 6)
        expected = self.read_data(self.dataPath+"/minimum_skew_outputs/output_", "int", 6)
        actual = []
        for input in inputs:
            actual.append(self.uut.skew_min(input[0]))

        self.assertEqual(expected, actual)

    def test_approximate_pattern_match(self):
        inputs = self.read_data(self.dataPath+"/approximate_pattern_match_inputs/input_", "str", 8)
        expected = self.read_data(self.dataPath+"/approximate_pattern_match_outputs/output_", "int", 8)
        actual = []
        for input in inputs:
            actual.append(self.uut.approximate_pattern_match(input[1], input[0], int(input[2])))

        self.assertEqual(expected, actual)

    def test_approximate_pattern_count(self):
        inputs = self.read_data(self.dataPath+"/approximate_pattern_count_inputs/input_", "str", 2)
        expected = self.read_data(self.dataPath+"/approximate_pattern_count_outputs/output_", "int", 2)
        actual = []
        for input in inputs:
            actual.append([self.uut.approximate_pattern_count(input[0], input[1], int(input[2]))])

        self.assertEqual(expected, actual)


if __name__ == '__main__':
    unittest.main()
