
from collections import deque
from typing import List, OrderedDict


class AssemblyAlgorithm:

    def __init__(self) -> None:
        pass

    def composition(self, k: int, text: str) -> List[str]:
        """ Composition

        Given a string Text, its k-mer composition Compositionk(Text) is the
        collection of all k-mer substrings of Text (including repeated k-mers).

        For example,
        Composition(3, TATGGGGTGC) = {ATG, GGG, GGG, GGT, GTG, TAT, TGC, TGG}.
        """

        res = []
        for i in range(len(text)-k+1):
            res.append(text[i:i+k])

        return res

    def genome_path(self, kmers: List[str]) -> str:
        """ Genome Path

        Reconstruct a string from its genome path.
        Input:  A sequence path of k-mers Pattern_1, ... ,Pattern_n such that
                the last k - 1 symbols of Pattern_i are equal to the first k-1
                symbols of Pattern_i+1 for 1 <= i <= n-1.
        Output: A string Text of length k+n-1 such that the i-th k-mer in Text
                is equal to Pattern_i (for 1 <= i <= n).
        """
        k = len(kmers[0])

        genome = kmers[0]
        for i in range(1, len(kmers)):
            print(kmers[i], kmers[i][k-1])
            genome += kmers[i][k-1]

        return genome

    def prefix(self, s: str, k: int) -> str:
        return s[:k]

    def suffix(self, s: str, k: int) -> str:
        return s[-k:]

    def overlap(self, kmers: List[str]) -> List[List[str]]:
        """ Overlap Graph

        Input:  A collection Patterns of k-mers.
        Output: The overlap graph Overlap(Patterns), in the form of an adjacency
                list.

        Note: Does not account for repeated elements in Patterns.
        """
        graph = []
        for kmer in kmers:
            graph.append([kmer])

        k = len(kmers[0]) - 1
        n = len(kmers)
        for i in range(n):
            for j in range(n):
                if j == i:
                    continue

                if self.suffix(kmers[i], k) == self.prefix(kmers[j], k):
                    if kmers[j] not in graph[i]:
                        graph[i].append(kmers[j])

        return graph

    def print_adjacency_graph(self, graph: List[List[str]], discardNoEdge=True):

        for adj in graph:
            if discardNoEdge and len(adj) < 2:
                continue
            elif not discardNoEdge and len(adj) < 2:
                print(adj[0])
                continue

            s = adj[0] + "->"
            for i in range(1, len(adj)):
                s += adj[i] + ","
            print(s[:-1])

    def path_graph(self, k: int, text: str) -> List[str]:
        """ Path Graph

        PathGraph(k, Text) is the path consisting of |Text| - k + 1 edges, where
        the i-th edge of this path is labeled by the i-th k-mer in Text and the
        i-th node of the path is labeled by the i-th (k - 1)-mer in Text
        """
        path = []
        kmers = self.composition(k, text)
        for kmer in kmers:
            path.append(self.prefix(kmer, k-1))
        path.append(self.suffix(kmers[-1], k-1))

        return path

    def deBruijn_graph(self, k: int, text: str) -> List[List[str]]:
        """ de Burijn Graph

        Given a genome Text, PathGraph(k, Text) is the path consisting of
        |Text| - k + 1 edges, where the i-th edge of this path is labeled by
        the i-th k-mer in Text and the i-th node of the path is labeled by the
        i-th (k - 1)-mer in Text. The de Bruijn graph DeBruijn(k, Text) is
        formed by gluing identically labeled nodes in PathGraphk(Text).
        """
        path = self.path_graph(k, text)

        map = OrderedDict()
        for kmer in sorted(path):
            map[kmer] = []

        for i in range(len(path)-1):
            map[path[i]].append(path[i+1])

        graph = []
        for key, val in map.items():
            graph.append([key] + val)

        return graph

    def deBruijn_graph_pattern(self, kmers: List[str]) -> List[List[str]]:
        """ de Bruijn Graph

        Construct de Bruijn graphs without gluing. Given a collection of k-mers
        Patterns, the nodes of DeBruijn(k, Patterns) are simply all unique
        (k−1)-mers occurring as a prefix or suffix in Patterns.

        For every k-mer in Patterns, we connect its prefix node to its suffix
        node by a directed edge in order to produce DeBruijn(Patterns).
        """
        k = len(kmers[0]) - 1

        map = OrderedDict()
        for kmer in sorted(kmers):
            p = self.prefix(kmer, k)
            s = self.suffix(kmer, k)
            if p in map:
                map[p].append(s)
            else:
                map[p] = [s]

        graph = []
        for key, val in map.items():
            graph.append([key] + sorted(val))

        return graph


    def make_adjacency_graph(self, input: List[str]) -> dict:

        graph = {}
        for adjacency in sorted(input):
            tokens = adjacency.split("->")
            node = tokens[0].strip()
            connectedTo = sorted([x.strip() for x in tokens[1].split(",")])
            graph[node] = connectedTo
        return graph





    def eulerian_path(self, graph: dict) -> List:

        edges = []
        for node, adj_nodes in graph.items():
            for node_j in adj_nodes:
                edges.append((node, node_j))

        print(edges)
        visited = []
        cycles = []
        while edges:

            edge = edges.pop(0)
            if edge in visited:
                continue

            path = [edge[0]]
            cycle = []
            while path:
                node = path.pop()
                if cycle and node == cycle[0]:
                    break

                cycle.append(node)
                for adj_node in graph[node]:
                    if (node, adj_node) not in visited:
                        path.append(adj_node)
                        visited.append((node, adj_node))
                        break


            cycles.append(cycle)
            print(cycles)

        return self.combine_cycles(cycles)

    def combine_cycles(self, cycles: List[List]) -> List:

        path = cycles.pop()

        while cycles:

            sub_path = cycles.pop()
            intersect = (set(path) & set(sub_path)).pop()

            print(intersect, path, sub_path)

            combined = []
            for i in range(sub_path.index(intersect), len(sub_path)):
                combined.append(sub_path[i])
            for i in range(0, sub_path.index(intersect)):
                combined.append(sub_path[i])
            for i in range(path.index(intersect), len(path)):
                combined.append(path[i])
            for i in range(0, path.index(intersect)):
                combined.append(path[i])

            path = combined

        return path

uut = AssemblyAlgorithm()

input = [
    "1 -> 2",
    "2 -> 3",
    "3 -> 4",
    "4 -> 5,7",
    "5 -> 6",
    "6 -> 4",
    "7 -> 8",
    "8 -> 9",
    "9 -> 7"
]
graph = uut.make_adjacency_graph(input)


# print(uut.eulerian_path(graph))


l1 = [0,3,2,1]
l2 = [2,6,5,4]
l3 = [6,8,7,9]


print(uut.combine_cycles([l1, l2, l3]))
