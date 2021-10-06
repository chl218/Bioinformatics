import unittest
import numpy as np
from src.tba import TBA


class TestTBA(unittest.TestCase):

    def __init__(self, methodName: str = ...) -> None:
        super().__init__(methodName=methodName)

        self.uut = TBA()
        self.dataPath = "01_Finding_Hidden_Messages_In_DNA/data"

    def read_data(self, file_path: str, read_type: str, file_count: int) -> list:
        data = []
        for i in range(1, file_count):
            res = np.loadtxt(file_path+str(i)+".txt", dtype=read_type)
            if res.size == 1:
                data.append([res.tolist()])
            else:
                data.append(res.tolist())
        return data

    def test_pattern_count(self):

        inputs = self.read_data(
            self.dataPath+"/pattern_count_inputs/input_", "str", 9)
        expected = self.read_data(
            self.dataPath+"/pattern_count_outputs/output_", "int", 9)
        actual = []
        for input in inputs:
            actual.append([self.uut.pattern_count(input[0], input[1])])

        self.assertEqual(expected, actual)

    def test_frequent_words(self):

        inputs = self.read_data(
            self.dataPath+"/frequent_words_inputs/input_", "str", 6)
        expected = self.read_data(
            self.dataPath+"/frequent_words_outputs/output_", "str", 6)
        actual = []
        for input in inputs:
            actual.append(self.uut.better_frequent_words(
                input[0], int(input[1])))

        for a1, a2 in zip(expected, actual):
            self.assertCountEqual(a1, a2)


if __name__ == '__main__':
    unittest.main()
