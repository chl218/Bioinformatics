import numpy as np
from HammingDistance import HammingDistance

def ApproximatePatternMatching(Text, Pattern, d):
    """
    Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.
        Input: Strings Pattern and Text along with an integer d.
        Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.
     
    Sample Input:
        CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT
        ATTCTGGA
        3
    Sample Output:
        6 7 26 27

     """
    idx = []
    map = {}    # Use map to avoid redundent Hamming Distance computation
    for i in range(0, len(Text) - len(Pattern) + 1):
        q = Text[i:i+len(Pattern)]
        if q in map and map[q] <= d:
            idx.append(i)
        else:
            map[q] = HammingDistance(Pattern, q)
            if map[q] <= d:
                idx.append(i)
                
    return idx

# %%
data = np.loadtxt('./data/dataset_9_4.txt', dtype='str')
print('Text: {0}\n\nPattern: {1}\n\nApproximate Pattern Matching Indicies: {2}'. 
      format(data[1], data[0],ApproximatePatternMatching(data[1], data[0], int(data[2]))))