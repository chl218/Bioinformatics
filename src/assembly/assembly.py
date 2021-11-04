
from typing import List


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

        path = []
        kmers = self.composition(k, text)
        for kmer in kmers:
            path.append(self.prefix(kmer, k-1))
        path.append(self.suffix(kmers[-1], k-1))

        return path

uut = AssemblyAlgorithm()

# input = None
# with open("/home/chl218/Downloads/dataset_198_10.txt") as f:
#     input = f.read().splitlines()
# uut.print_adjacency_graph(uut.overlap(input))



# k = 4
# t = "AAGATTCTCTAAGA"
# kmers = uut.composition(k, t)
# print(t, kmers)

# uut.print_adjacency_graph(uut.overlap(kmers), False)


k = 3
t = "TAATGCCATGGGATGTT"

print(uut.path_graph(3, t))