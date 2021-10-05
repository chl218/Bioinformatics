import unittest
import numpy as np
from src.tba import tba

class TestTba(unittest.TestCase):

    def __init__(self, methodName: str = ...) -> None:
        super().__init__(methodName=methodName)

        self.uut = tba()
        self.dataPath = "01-Finding-Hidden-Messages-in-DNA/data"

    def test_pattern_count(self):

        fin = []
        expected = []

        for i in range(1,9):
           fin.append(np.loadtxt(self.dataPath+"/pattern-count-inputs/input_"+str(i)+".txt", dtype='str'))
           expected.append(np.loadtxt(self.dataPath+"/pattern-count-outputs/output_"+str(i)+".txt", dtype='int').item(0))

        actual = []
        for f in fin:
            actual.append(self.uut.pattern_count(f[0], f[1]))

        self.assertEqual(expected, actual)


if __name__ == '__main__':
    unittest.main()
