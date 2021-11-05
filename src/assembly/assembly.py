
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



    def deBruijn_graph_kmers(self, kmers: List[str]) -> List[List[str]]:

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

uut = AssemblyAlgorithm()

input = None
with open("/home/chl218/Downloads/dataset_200_8.txt") as f:
    input = f.read().splitlines()
# uut.print_adjacency_graph(uut.overlap(input))


kmers = [
    "GAGG",
    "CAGG",
    "GGGG",
    "GGGA",
    "CAGG",
    "AGGG",
    "GGAG"
]

uut.print_adjacency_graph(uut.deBruijn_graph_kmers(input))