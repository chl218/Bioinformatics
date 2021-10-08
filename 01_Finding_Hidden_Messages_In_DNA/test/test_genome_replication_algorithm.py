import unittest
import numpy as np
from src.genome_replication_algorithm import GenomeReplicationAlgorithm


class TestGenomeReplicationAlgorithm(unittest.TestCase):

    def __init__(self, methodName: str = ...) -> None:
        super().__init__(methodName=methodName)

        self.uut = GenomeReplicationAlgorithm()
        self.dataPath = "01_Finding_Hidden_Messages_In_DNA/data"

    def read_data(self, file_path: str, read_type: str, file_count: int) -> list:
        data = []
        for i in range(1, file_count+1):
            res = np.loadtxt(file_path+str(i)+".txt", dtype=read_type)
            if res.size == 1:
                data.append([res.tolist()])
            else:
                data.append(res.tolist())
        return data