import numpy as np

def Skew(genome):
    """
    We define Skewi(Genome) as the difference between the total number of occurrences of G and the total number of
    occurrences of C in the first i nucleotides of Genome.
    """
    
    difference = 0
    skewArray = [0]
    
    for c in genome:
        if c == 'G':
            difference += 1
            skewArray.append(difference)
        elif c == 'C':
            difference -= 1
            skewArray.append(difference)
        else:
            skewArray.append(difference)
            
    return skewArray
    

def MinimumSkew(genome):
    """
    Minimum Skew Problem: Find a position in a genome minimizing the skew.
        Input:  A DNA string Genome.
        Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).
        
    Sample Input:
        TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT
    Sample Output:
        11 24
    """
    skewArr = np.array(Skew(genome))
    return np.where(skewArr==skewArr.min())[0]


# %%
data = np.loadtxt('./data/dataset_7_6.txt', dtype='str', ndmin=1)
print('DNA: {0}\n\nMinimum Skew Index: {1}'.format(data[0], MinimumSkew(data[0])))
