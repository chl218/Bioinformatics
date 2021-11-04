from typing import Any, List
import unittest
import numpy as np
from src.assembly.assembly import AssemblyAlgorithm


class TestAssemblyAlgorithm(unittest.TestCase):

    def __init__(self, methodName: str = ...) -> None:
        super().__init__(methodName=methodName)

        self.uut = AssemblyAlgorithm()
        self.dataPath = "data/assembly"


    def read_data(self, file_path: str, read_type: str, file_count: int) -> List:

        data = []
        for i in range(file_count):
            with open(file_path+str(i)+".txt") as f:
                if read_type == "int":
                    line = f.read().splitlines()
                    data.append(list(map(int, line)))
                else:
                    data.append(f.read().splitlines())

        return data

    def input_to_graph(self, xs: List[str]) -> List[List[str]]:

        graph = []
        sorted(xs)
        for x in sorted(xs):
            nodes = x.split("->")

            node = nodes[0]
            connects = sorted(nodes[1].split(","))

            graph.append([node] + connects)

        return graph


    def compare_graphs(self, graph1: List[List[str]], graph2: List[List[str]]):
        for g1, g2 in zip(graph1, graph2):
            for adj_lst in g2:
                if len(adj_lst) == 1:
                    continue
                self.assertIn(adj_lst, g1)


    def test_composition(self):
        inputs = self.read_data(self.dataPath+"/composition_inputs/test", "str", 4)
        expected = self.read_data(self.dataPath+"/composition_outputs/test", "str", 4)
        actual = []
        for input in inputs:
            actual.append(self.uut.composition(int(input[0]), input[1]))

        self.assertEqual(expected, actual)

    def test_overlap_graph(self):
        inputs = self.read_data(self.dataPath+"/overlap_graph_inputs/test", "str", 1)
        expected = self.read_data(self.dataPath+"/overlap_graph_outputs/test", "str", 1)
        actual = []
        for input in inputs:
            actual.append(self.uut.overlap(input))



        # # self.assertEqual(expected, actual)

        print(inputs[0])
        print(expected[0])
        print(self.input_to_graph(expected[0]))
        print(actual[0])

        for g1, g2 in zip(expected, actual):
            g1 = self.input_to_graph(g1)
            for adj_lst in g2:
                if len(adj_lst) == 1:
                    continue
                self.assertIn(adj_lst, g1)

