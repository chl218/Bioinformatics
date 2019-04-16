import numpy as np
from Neighbors import Neighbors

def FrequentWordsWithMismatches(Text, k, d):
    """
    Frequent Words with Mismatches Problem: Find the most frequent k-mers with mismatches in a string.
        Input: A string Text as well as integers k and d. (You may assume k ≤ 12 and d ≤ 3.)
        Output: All most frequent k-mers with up to d mismatches in Text.
     
     
    Sample Input:
        ACGTTGCATGTCGCATGATGCATGAGAGCT 
        4
        1
    Sample Output:
        ATGC ATGT GATG
    """
    frequentPatterns = set()
    frequencyArray = ApproximateFrequencies(Text, k, d)
    maxCount = max(frequencyArray.values())
    
    for key in frequencyArray:
        if frequencyArray[key] == maxCount:
            frequentPatterns.add(key)
            
    return frequentPatterns


def ApproximateFrequencies(Text, k, d):
    """
    Count the k-mers in text, with d Hamming Distance
    """
    frequencyArray = {}
    
    for i in range(0, len(Text) - k + 1):
        p = Text[i:i+k]
        neighbors = Neighbors(p, d)
        for neighbor in neighbors:
            if neighbor in frequencyArray:
                frequencyArray[neighbor] += 1
            else:
                frequencyArray[neighbor] = 1
        
    return frequencyArray

# %%
data = np.loadtxt('./data/dataset_9_7.txt', dtype='str')
print('Text: {0}\n\nk: {1}\n\nd: {2}\n\nMost Frequent Words With Mismatches: {3}'
      .format(data[0], data[1], data[2], FrequentWordsWithMismatches(data[0], int(data[1]), int(data[2]))))