import unittest
import numpy as np
from src.assembly.assembly import AssemblyAlgorithm


class TestAssemblyAlgorithm(unittest.TestCase):

    def __init__(self, methodName: str = ...) -> None:
        super().__init__(methodName=methodName)

        self.uut = AssemblyAlgorithm()
        self.dataPath = "data/assembly"


    def read_data(self, file_path: str, read_type: str, file_count: int) -> list:

        data = []
        for i in range(file_count):
            with open(file_path+str(i)+".txt") as f:
                if read_type == "int":
                    line = f.read().splitlines()
                    data.append(list(map(int, line)))
                else:
                    data.append(f.read().splitlines())

        return data

    def test_composition(self):
        inputs = self.read_data(self.dataPath+"/composition_inputs/test", "str", 4)
        expected = self.read_data(self.dataPath+"/composition_outputs/test", "str", 4)
        actual = []
        for input in inputs:
            actual.append(self.uut.composition(int(input[0]), input[1]))

        self.assertEqual(expected, actual)
