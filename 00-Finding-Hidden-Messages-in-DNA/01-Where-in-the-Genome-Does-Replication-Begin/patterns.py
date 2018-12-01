# %%
import numpy as np


# %%

# Search the amount of pattern matched in text
def PatternCount(Text, Pattern):

    textLen    = len(Text)
    patternLen = len(Pattern)
    
    count = 0
    for i in range(0, textLen - patternLen):
        if Text[i:i+patternLen] == Pattern:
            count += 1
    
    return count

# Find the most frequent k-mer in text
def FrequentWords(Text, k):
    textLen = len(Text)
    
    frequentPatterns = set()
    count = [0] * textLen
    
    # Search for k-mer patterns
    for i in range(0, textLen - k):
        pattern = Text[i:i+k]
        count[i] = PatternCount(Text, pattern)
    maxCount = max(count)
    
    # Search for most frequent k-mers
    for i in range(0, textLen - k):
        if count[i] == maxCount:
            frequentPatterns.add(Text[i:i+k])
            
    return frequentPatterns

# %%    
genome = np.loadtxt('dataset_2_7.txt', dtype='str');
print(PatternCount(genome[0], genome[1]))

# %%
genome = np.loadtxt('dataset_2_10.txt', dtype='str');
res = FrequentWords(genome[0], int(genome[1]))
for seq in res:
    print(seq, ' ', end='')