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
    
    # Search for most frequent k-mers
    for i in range(0, textLen - k + 1):
        if count[i] == maxCount:
            frequentPatterns.add(Text[i:i+k])
            
    return frequentPatterns

# Return the reversed complement of DNA pattern
def ReverseComplement(strand):
    res = []
    
    for p in strand:
        if p == 'T':
            res.append('A')
        elif p == 'A':
            res.append('T')
        elif p == 'G':
            res.append('C')
        elif p == 'C':
            res.append('G')
            
    return ''.join(res[::-1])

# Return the index in the genome that matches to the given pattern
def PatternMatching(pattern, genome):
    # regular expression with lookahead e.g. re.finditer((?=pattern), text)
    return [iter.start() for iter in re.finditer('(?='+pattern+')', genome)]
    
# Find all k-mer in genome
def FindAllKmers(genome, k):
    textLen = len(genome)
    
    frequentPatterns = set()
    count = [0] * textLen
    
    # Search for k-mer patterns
    for i in range(0, textLen - k + 1):
        pattern = genome[i:i+k]
        count[i] = PatternCount(genome, pattern)
    # Search for most frequent k-mers
    for i in range(0, textLen - k + 1):
        if count[i] > 1:
            frequentPatterns.add(genome[i:i+k])
            
    return frequentPatterns

# Find all k-mer appearence >= t in genome of length L
def ClumpFinding(genome, k, L, t):
    clumps = []
    
    # print(k, L, t)
    kmers = FindAllKmers(genome, k)
    # print(kmers)
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
            #print('i =',i, count, displacement)
            for j in range(i+1, len(idx)):
                #print('j =',j, count, displacement, idx[j])
                if(displacement - idx[j] + k > L):
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



# %% 

def PatternToNumber(pattern):
    num = 0
    k   = 0
    for p in pattern[::-1]:
        if   p == 'A':  num += 0 * 4**k
        elif p == 'C':  num += 1 * 4**k
        elif p == 'G':  num += 2 * 4**k
        elif p == 'T':  num += 3 * 4**k
        k += 1
        
    return num

def NumberToPattern(number, k):
    pattern = ""
    
    for i in range(0, k):
        r = number % 4
        number = int(number/4)

        if   r == 0: pattern = 'A' + pattern
        elif r == 1: pattern = 'C' + pattern
        elif r == 2: pattern = 'G' + pattern
        elif r == 3: pattern = 'T' + pattern    

    return pattern

# Return the frequence occurance of k-mer in genome
def ComputingFrequencies(genome, k):
    frequencyArray = {}
    
    for i in range(0, len(genome) - (k-1)):
        key = genome[i:i+k]
        if key in frequencyArray:
            frequencyArray[key] += 1
        else:
            frequencyArray[key] = 1
            
    return frequencyArray
    
# Return the maximum frequence occurance of k-mers in genome
def FasterFrequentWords(genome, k):
    frequentPatterns = set()
    frequencyArray = ComputingFrequencies(genome, k)
    maxCount = max(frequencyArray.values())
    
    for key in frequencyArray:
        if frequencyArray[key] == maxCount:
            frequentPatterns.add(NumberToPattern(key, k))
            
    return frequentPatterns
    
# Find all k-mer appearence >= t in genome of length L
def BetterClumpFinding(genome, k, L, t):
    clumps = set()
    # print(k, t, L)
    text = genome[0:L]
    frequencyMap = ComputingFrequencies(text, k)
    
    for key in frequencyMap:
        if frequencyMap[key] >= t:
            clumps.add(key)
            
    for i in range(1, len(genome) - L):
        prevPattern = genome[i-1:i-1+k]
        if prevPattern in frequencyMap:
            frequencyMap[prevPattern] -= 1
            
        nextPattern = genome[i+L-k:i+L]
        if nextPattern in frequencyMap:
            frequencyMap[nextPattern] += 1
        else:
            frequencyMap[nextPattern] = 1
            
        # print(i-1,i-1+k,prevPattern, i+L-k,i+L,nextPattern)
        if frequencyMap[nextPattern] >= t:
            clumps.add(nextPattern)
            
    return clumps


# %%    
genome = np.loadtxt('dataset_2_7.txt', dtype='str');
print('Pattern Count: ' + str(PatternCount(genome[0], genome[1])))

# %%
genome = np.loadtxt('dataset_2_10.txt', dtype='str');
print('Frequent Words: ', end='')
print(FrequentWords(genome[0], int(genome[1])))

    
# %%

g = np.loadtxt('dataset_3_2.txt', dtype='str', ndmin=1)
print('Reverse Complement: ', end='')
print(ReverseComplement(g[0]))

# %%

g = np.loadtxt('dataset_3_5.txt', dtype='str')
print('Pattern Matching: ', end='')
for i in PatternMatching(g[0], g[1]):
    print(i, '', end='')
    
# %%
    
vibrio_cholerae = np.loadtxt('Vibrio_cholerae.txt', dtype='str', ndmin=1)
for i in PatternMatching('CTTGATCAT', vibrio_cholerae[0]):
    print(i, '', end='')
    

# %%
data = np.loadtxt('dataset_4_5.txt', dtype='str', delimiter='\n ')
param = [int(s) for s in data[1].split(' ')]
print('Clump Finding: ', end='')
for val in ClumpFinding(data[0], param[0], param[1], param[2]):
    print(val, '', end='')


# %%
    
ecoli = np.loadtxt('E_coli.txt',  dtype='str', ndmin=1)
print(len(BetterClumpFinding(ecoli[0], 9, 500, 3)))
    

# %%

import numpy as np

def Skew(genome):
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

    
def MinSkew(genome):
    skewArr = np.array(Skew(genome))
    return np.where(skewArr==skewArr.min())[0]


def HammingDistance(p, q):
    if(len(p) != len(q)):
        return float('inf')
    
    count = 0;
    for c1, c2 in zip(p, q):
        if c1 != c2:
            count += 1
            
    return count


def ApproxPatternMatching(pattern, genome, distance):
    idx = []
    map = {}
    cnt = 0
    for i in range(0, len(genome) - len(pattern) + 1):
        q = genome[i:i+len(pattern)]
        if q in map and map[q] <= distance:
            idx.append(i)
            cnt += 1
        else:
            map[q] = HammingDistance(pattern, q)
            if map[q] <= distance:
                idx.append(i)
                cnt += 1
                
    return idx, cnt



# %%
    
genome = np.loadtxt('dataset_7_6.txt', dtype='str', ndmin=1)
print(MinSkew(genome[0]))

# %%

genomes = np.loadtxt('dataset_9_3.txt', dtype='str')
HammingDistance(genomes[0], genomes[1])

# %%
ApproxPatternMatching('ATTCTGGA','CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT',3)


data = np.loadtxt('dataset_9_4.txt', dtype='str')
for i in ApproxPatternMatching(data[0], data[1], int(data[2])):
    print(i, end=' ')

# %%
data = np.loadtxt('dataset_9_6.txt', dtype='str')
_, cnt = ApproxPatternMatching(data[0],data[1],int(data[2]))
cnt