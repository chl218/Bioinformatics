from typing import List


class GenomeReplicationAlgorithm:

    def __init__(self) -> None:
        pass

    def skew_GC(self, genome: str) -> List[int]:
        """Skew Diagram

        Skew diagram between G C nucleotides
        """

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

    def skew_min(self, genome: str) -> List[int]:
        """Minimum Skew

        Find the position(s) in a genome where the skew diagram attains a minimum.
        """

        skew = self.skew_GC(genome)
        mmin = min(skew)
        return [idx for idx, val in enumerate(skew) if val == mmin]


uut = GenomeReplicationAlgorithm()
