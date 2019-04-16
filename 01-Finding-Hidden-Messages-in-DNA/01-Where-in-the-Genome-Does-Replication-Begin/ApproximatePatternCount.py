import numpy as np
from HammingDistance import HammingDistance

def ApproximatePatternCount(Text, Pattern, d):
    """
    ApproximatePatternCount Problem:
        Input:  Strings Pattern and Text, and an integer d
        Output: The number of times Pattern appears in Text with at most d mismatches
  
    Sample Input:
        GAGG
        TTTAGAGCCTTCAGAGG
        2
    Sample Output:
        4
    """
    map = {}
    count = 0
    for i in range(0, len(Text) - len(Pattern) + 1):
        q = Text[i:i+len(Pattern)]
        if q in map and map[q] <= d:
            count += 1
        else:
            map[q] = HammingDistance(Pattern, q)
            if map[q] <= d:
                count += 1
    return count


# %%
data = np.loadtxt('./data/dataset_9_6.txt', dtype='str')
print('Pattern: {0}\n\nText: {1}\n\nDistance: {2}\n\nApproximate Pattern Count: {3}'.
      format(data[0], data[1], data[2], ApproximatePatternCount(data[1],data[0],int(data[2]))))

