import numpy as np
from HammingDistance import HammingDistance

def Neighbors(Pattern, d):
    """
    Neighbors Problem:
        Input: A string Pattern and an integer d.
        Output: The collection of strings Neighbors(Pattern, d).
     
        
    Sample Input:
        ACG
        1
    Sample Output:
        AAG ACG CCG AGG ACT GCG ATG ACC TCG ACA
    """

    if d == 0:
        return [Pattern]
    if len(Pattern) == 1:
        return ['A','C','G','T']
    
    neighbors = set()
    for substr in Neighbors(Pattern[1:], d):
        if HammingDistance(Pattern[1:], substr)  < d:
            for nucleotide  in ['A','C','G','T']:
                neighbors.add(nucleotide + substr)
        else:
            neighbors.add(Pattern[0] + substr)
            
    return neighbors

# %%
    
data = np.loadtxt('./data/dataset_3014_4.txt', dtype='str')
print('Pattern: {0}\n\nHamming Distance: {1}\n\nNeighbors: {2}'
      .format(data[0], data[1], Neighbors(data[0], int(data[1]))))

    