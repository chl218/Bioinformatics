import numpy as np
from ApproximatePatternCount import ApproximatePatternCount
from FrequentWordsWithMismatches import ApproximateFrequencies
from ReverseComplement import ReverseComplement

def FrequentWordsWithMismatchesAndReverseComplements(Text, k, d):
    """
    Frequent Words with Mismatches and Reverse Complements Problem: Find the most frequent k-mers (with mismatches and
    reverse complements) in a DNA string.
        Input:  A DNA string Text as well as integers k and d.
        Output: All k-mers Pattern maximizing the sum Countd(Text, Pattern) + Countd(Text, Pattern) over all possible
                k-mers.
    
    
    Sample Input:
        ACGTTGCATGTCGCATGATGCATGAGAGCT
        4
        1
    Sample Output:
        ACAT ATGT
    """
    freqPattern = set()
    
    freqArr = ApproximateFrequencies(Text, k, d)
    
    sumFreqArr = {}
    for p in freqArr:
        pp = ReverseComplement(p)
        
        if (p, pp) not in sumFreqArr and (pp, p) not in sumFreqArr:
            count = ApproximatePatternCount(Text, pp, d)
            sumFreqArr[(p, pp)] = freqArr[p] + count

    maxCount = max(sumFreqArr.values())
    
    for p in sumFreqArr:
        if sumFreqArr[p] == maxCount:
            freqPattern.add(p[0])
            freqPattern.add(p[1])
        
    return freqPattern


# %%
data = np.loadtxt('./data/dataset_9_8.txt', dtype='str')
print('Text: {0}\n\nk: {1}\n\nd: {2}\n\nMost Frequent Words With Mismatches And Reverse Complements: {3}'
      .format(data[0], data[1], data[2], FrequentWordsWithMismatchesAndReverseComplements(data[0], int(data[1]), int(data[2]))))