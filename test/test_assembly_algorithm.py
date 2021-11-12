from typing import List
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
        for x in sorted(xs):
            tokens = x.split("->")
            node = tokens[0].strip()
            connectedTo = sorted([x.strip() for x in tokens[1].split(",")])
            graph.append([node] + connectedTo)

        return graph

    def compare_graphs(self, expected_graph: List[List[str]], actual_graph: List[List[str]]):
        for adj_lst in actual_graph:
            if len(adj_lst) == 1:
                continue
            adj_lst[1:] = sorted(adj_lst[1:])
            self.assertIn(adj_lst, expected_graph)

    def test_composition(self):
        inputs = self.read_data(self.dataPath+"/composition_inputs/test", "str", 5)
        expected = self.read_data(self.dataPath+"/composition_outputs/test", "str", 5)
        actual = []
        for input in inputs:
            actual.append(self.uut.composition(int(input[0]), input[1]))

        self.assertEqual(expected, actual)

    def test_overlap_graph(self):
        inputs = self.read_data(self.dataPath+"/overlap_graph_inputs/test", "str", 7)
        expected = self.read_data(self.dataPath+"/overlap_graph_outputs/test", "str", 7)
        actual = []
        for input in inputs:
            actual.append(self.uut.overlap(input))

        for g1, g2 in zip(expected, actual):
            self.compare_graphs(self.input_to_graph(g1), g2)

    def test_deBruijn_graph(self):
        inputs = self.read_data(self.dataPath+"/deburijn_graph_inputs/test", "str", 5)
        expected = self.read_data(self.dataPath+"/deburijn_graph_outputs/test", "str", 5)
        actual = []
        for input in inputs:
            actual.append(self.uut.deBruijn_graph(int(input[0]), input[1]))

        for g1, g2 in zip(expected, actual):
            self.compare_graphs(self.input_to_graph(g1), g2)

    def test_deBruijn_graph_pattern(self):
        inputs = self.read_data(self.dataPath+"/deburijn_graph_patterns_inputs/test", "str", 5)
        expected = self.read_data(self.dataPath+"/deburijn_graph_patterns_outputs/test", "str", 5)
        actual = []
        for input in inputs:
            actual.append(self.uut.deBruijn_graph_pattern(input))

        for g1, g2 in zip(expected, actual):
            self.compare_graphs(self.input_to_graph(g1), g2)

    def test_eulerian_cycle(self):
        expected = self.read_data(self.dataPath+"/eulerian_cycle_outputs/test", "str", 7)

        inputs = []
        for input in self.read_data(self.dataPath+"/eulerian_cycle_inputs/test", "str", 7):
            inputs.append(self.uut.make_adjacency_graph(input))

        actual = []
        for input in inputs:
            actual.append(self.uut.eulerian_cycle(input))

        for cycle, graph in zip(actual, inputs):
            edges = []
            for src, dsts in graph.items():
                for dst in dsts:
                    edges.append((src, dst))

            visited_edges = 0
            for i in range(len(cycle)-1):
                self.assertIn((cycle[i], cycle[i+1]), edges)
                visited_edges += 1

            self.assertEqual(len(edges), visited_edges)

    def test_eulerian_path(self):
        expected = self.read_data(self.dataPath+"/eulerian_path_outputs/test", "str", 6)

        inputs = []
        for input in self.read_data(self.dataPath+"/eulerian_path_inputs/test", "str", 6):
            inputs.append(self.uut.make_adjacency_graph(input))

        actual = []
        for input in inputs:
            actual.append(self.uut.eulerian_path(input))

        for cycle, graph in zip(actual, inputs):
            edges = []
            for src, dsts in graph.items():
                for dst in dsts:
                    edges.append((src, dst))

            visited_edges = 0
            for i in range(len(cycle)-1):
                self.assertIn((cycle[i], cycle[i+1]), edges)
                visited_edges += 1

            self.assertEqual(len(edges), visited_edges+1)


if __name__ == '__main__':
    unittest.main()
