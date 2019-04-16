# %%
def HammingDistance(p, q):
    """
    Hamming Distance Problem: Compute the Hamming distance between two strings.
        Input: Two strings of equal length.
        Output: The Hamming distance between these strings.
        
    Sample Input:
        GGGCCGTTGGT
        GGACCGTTGAC
    Sample Output:
        3
    """
    count = 0;
    for c1, c2 in zip(p, q):
        if c1 != c2:
            count += 1
    return count

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


# Find all k-mer in genome
def FindAllKmers(dna, k):
    textLen = len(dna)
      
    kmers = set()
    # Search for k-mer patterns
    for i in range(0, textLen - k + 1):
        kmers.add(dna[i:i+k])
        
    return kmers

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
def MotifEnumeration(dna, k, d):
    """
    Motif Enumeration
        Input: Integers k and d, followed by a collection of string DNA
        Output: All (k, d)-motifs in DNA
    """
    
    patterns = set()
    
    for _dna in dna:
        kmers = FindAllKmers(_dna, k)
        for kmer in kmers:
            _patterns = Neighbors(kmer, d) # find all pattern with d mismatch
            for pattern in _patterns:
                
                count = 0;
                for __dna in dna:
                    idx = ApproximatePatternMatching(__dna, pattern, d)
                    if idx.count > 0:
                        count += 1
                        
                if count == len(dna):
                    patterns.add(pattern)
    
    return patterns
# %%
        