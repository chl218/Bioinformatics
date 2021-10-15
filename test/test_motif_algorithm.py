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

if __name__ == '__main__':
    unittest.main()
