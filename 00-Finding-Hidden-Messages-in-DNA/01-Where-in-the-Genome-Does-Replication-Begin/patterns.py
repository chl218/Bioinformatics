# %%
import numpy as np
import re

# %%

# Search the amount of pattern matched in text
def PatternCount(Text, Pattern):

    textLen    = len(Text)
    patternLen = len(Pattern)
    
    count = 0
    for i in range(0, textLen - patternLen + 1):
        if Text[i:i+patternLen] == Pattern:
            count += 1
    
    return count

# Find the most frequent k-mer in text
def FrequentWords(Text, k):
    textLen = len(Text)
    
    frequentPatterns = set()
    count = [0] * textLen
    
    # Search for k-mer patterns
    for i in range(0, textLen - k + 1):
        pattern = Text[i:i+k]
        count[i] = PatternCount(Text, pattern)
    maxCount = max(count)
    
    print(count)
    # Search for most frequent k-mers
    for i in range(0, textLen - k + 1):
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

def PatternMatching(pattern, genome):
    # regular expression with lookahead e.g. re.finditer((?=pattern), text)
    return [iter.start() for iter in re.finditer('(?='+pattern+')', genome)]
    





# %%

# Find all k k-mer in text
def FindAllKmers(Text, k):
    textLen = len(Text)
    
    frequentPatterns = set()
    count = [0] * textLen
    
    # Search for k-mer patterns
    for i in range(0, textLen - k + 1):
        pattern = Text[i:i+k]
        count[i] = PatternCount(Text, pattern)
    # Search for most frequent k-mers
    for i in range(0, textLen - k + 1):
        if count[i] > 1:
            frequentPatterns.add(Text[i:i+k])
            
    return frequentPatterns

def ClumpFinding(genome, k, L, t):
    clumps = []
    
    print(k, L, t)
    
    kmers = FindAllKmers(genome, k)
    print(kmers)
    for kmer in kmers:
        idx = PatternMatching(kmer, genome)
    
        if(len(idx) < t):
            continue
    
        idx = idx[::-1]
        # print(kmer, idx)    
        
        for i in range(0, len(idx)):
            displacement = idx[i]
            count = 1
            flag = False
            # print('i =',i, count, displacement)

            for j in range(i+1, len(idx)):
                # print('j =',j, count, displacement, idx[j])

                if(displacement - idx[j] > L):
                    break
                else:
                    count += 1

                if(count == t):
                    clumps.append(kmer)
                    flag = True
                    break
                    
            if flag:
                break
        
    return clumps    

#data = np.loadtxt('dataset_4_5.txt', dtype='str', delimiter='\n ')
#param = [int(s) for s in data[1].split(' ')]
#print(ClumpFinding(data[0], param[0], param[1], param[2]))

print(ClumpFinding('CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA', 5, 50, 4))

print(ClumpFinding('AAAACGTCGAAAAA', 2, 4, 2))

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
for i in PatternMatching(g[0], g[1]):
    print(i, '', end='')
    
# %%
    
vibrio_cholerae = np.loadtxt('Vibrio_cholerae.txt', dtype='str', ndmin=1)
for i in PatternMatching('CTTGATCAT', vibrio_cholerae[0]):
    print(i, '', end='')
    


# %%
