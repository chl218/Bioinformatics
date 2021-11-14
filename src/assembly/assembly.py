
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

    def deBruijn_graph_pattern(self, kmers: List[str]) -> dict:
        """ de Bruijn Graph

        Construct de Bruijn graphs without gluing. Given a collection of k-mers
        Patterns, the nodes of DeBruijn(k, Patterns) are simply all unique
        (k−1)-mers occurring as a prefix or suffix in Patterns.

        For every k-mer in Patterns, we connect its prefix node to its suffix
        node by a directed edge in order to produce DeBruijn(Patterns).
        """
        k = len(kmers[0]) - 1

        graph = OrderedDict()
        for kmer in sorted(kmers):
            p = self.prefix(kmer, k)
            s = self.suffix(kmer, k)
            if p in graph:
                graph[p].append(s)
            else:
                graph[p] = [s]

        return graph

    def make_adjacency_graph(self, input: List[str]) -> dict:

        graph = {}
        for adjacency in sorted(input):
            tokens = adjacency.split("->")
            node = tokens[0].strip()
            connectedTo = sorted([x.strip() for x in tokens[1].split(",")])
            graph[node] = connectedTo
        return graph

    def eulerian_cycle(self, graph: dict) -> List:
        """ Eulerian Cycle

        Input: The adjacency list of an Eulerian directed graph.
        Output: An Eulerian cycle in this graph.
        """
        edges = []
        for node, adj_nodes in graph.items():
            for node_j in adj_nodes:
                edges.append((node, node_j))

        visited = []
        cycles = []

        # while there are unexplored edges in Graph
        while edges:

            #  select a node newStart in Cycle with still unexplored edges
            edge = edges.pop(0)
            if edge in visited:
                continue

            # form Cycle’ by traversing Cycle (starting at newStart) and then randomly walking
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

            # Cycle <- Cycle’
            cycles.append(cycle)

        return self.combine_cycles(cycles)

    def combine_cycles(self, cycles: List[List]) -> List:
        """ Combine Cycles

        Combine cycles in an Eulerian Graph, e.g. balanced graph where indegree
        and outdegree of a node v is the same for all nodes in the graph
        """

        path = cycles.pop()

        while cycles:

            sub_path = cycles.pop()
            intersect = (set(path) & set(sub_path))
            if not intersect:
                cycles.insert(0, sub_path)
                continue

            # Rotate the intersect node_v to front as the starting path, then
            # combine the two paths beginning and ending at node_v
            intersect = intersect.pop()
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

        return path + [path[0]]

    def eulerian_path(self, graph: dict) -> List:

        degrees = {}
        for node, adj_nodes in graph.items():
            if not node in degrees:
                degrees[node] = 0

            for node_j in adj_nodes:
                degrees[node] += 1
                if not node_j in degrees:
                    degrees[node_j] = -1
                else:
                    degrees[node_j] -= 1

        # Assume input graph is correct, there will only be two unbalanced nodes
        # since it is a dna string e.g. beginning/end of dna string
        src_node = None
        dst_node = None
        for node, degrees in degrees.items():
            if degrees > 0:
                dst_node = node
            if degrees < 0:
                src_node = node

        if not src_node in graph:
            graph[src_node] = [dst_node]
        else:
            graph[src_node].append(dst_node)


        # break the edge, rotate the path
        eulerian_cycle = self.eulerian_cycle(graph)
        src_idx = None
        dst_idx = None
        for i in range(len(eulerian_cycle)-1):
            if eulerian_cycle[i] == src_node and eulerian_cycle[i+1] == dst_node:
                src_idx = i
                dst_idx = i+1

        return eulerian_cycle[dst_idx:-1] + eulerian_cycle[:src_idx+1]

    def string_reconstruction(self, pattern: List[str]) -> str:

        graph = self.deBruijn_graph_pattern(pattern)
        path = self.eulerian_path(graph)
        return self.genome_path(path)


uut = AssemblyAlgorithm()


# k = 4
# pattern = [
#     "CTTA",
#     "ACCA",
#     "TACC",
#     "GGCT",
#     "GCTT",
#     "TTAC",
# ]

# print(uut.string_reconstruction(pattern))


# graph = [
#     "0 -> 2",
#     "1 -> 3",
#     "2 -> 1",
#     "3 -> 0,4",
#     "6 -> 3,7",
#     "7 -> 8",
#     "8 -> 9",
#     "9 -> 6"
# ]
# path = uut.eulerian_path(uut.make_adjacency_graph(graph))

# print('->'.join(path))

input = None
with open("/home/chl218/Downloads/dataset_203_7.txt") as f:
    input = f.read().splitlines()


print(uut.string_reconstruction(input[1:]))

# graph = uut.make_adjacency_graph(input)
# print('->'.join(uut.eulerian_path(graph)))


