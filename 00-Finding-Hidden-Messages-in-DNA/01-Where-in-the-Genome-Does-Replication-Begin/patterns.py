# %%
import numpy as np
import re

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

def ReverseComplement(s):
    res = []
    
    for p in s:
        if p == 'T':
            res.append('A')
        elif p == 'A':
            res.append('T')
        elif p == 'G':
            res.append('C')
        elif p == 'C':
            res.append('G')
            
    return ''.join(res[::-1])

def PatternIndex(pattern, genome):
    # regular expression with lookahead e.g. re.finditer((?=pattern), text)
    return [iter.start() for iter in re.finditer('(?='+pattern+')', genome)]
    
    
# %%    
genome = np.loadtxt('dataset_2_7.txt', dtype='str');
print(PatternCount(genome[0], genome[1]))

# %%
genome = np.loadtxt('dataset_2_10.txt', dtype='str');
res = FrequentWords(genome[0], int(genome[1]))
for seq in res:
    print(seq, ' ', end='')
    
# %%

g = np.loadtxt('dataset_3_2.txt', dtype='str', ndmin=1)
print(ReverseComplement(g[0]))

# %%

g = np.loadtxt('dataset_3_5.txt', dtype='str')
for i in PatternIndex(g[0], g[1]):
    print(i, '', end='')
    
# %%
    
vibrio_cholerae = np.loadtxt('Vibrio_cholerae.txt', dtype='str', ndmin=1)
for i in PatternIndex('CTTGATCAT', vibrio_cholerae[0]):
    print(i, '', end='')
    
