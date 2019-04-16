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