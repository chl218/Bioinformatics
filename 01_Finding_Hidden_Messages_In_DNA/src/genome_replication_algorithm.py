from typing import List


class GenomeReplicationAlgorithm:

    def __init__(self) -> None:
        pass

    def skew_GC(self, genome: str) -> List[int]:

        diff = 0
        skew = [0]

        n = len(genome)
        for i in range(0, n):
            if genome[i] == 'G':
                diff += 1
                skew.append(diff)
            elif genome[i] == 'C':
                diff -= 1
                skew.append(diff)
            else:
                skew.append(diff)

        return skew




uut = GenomeReplicationAlgorithm()
g = "GAGCCACCGCGATA"

print(' '.join(map(str, uut.skew_GC(g))))