# %%
import numpy as np

# %%
genome = np.loadtxt('dataset_2_7.txt', dtype='str');

# %%

def PatternCount(Text, Pattern):

    textLen    = len(Text)
    patternLen = len(Pattern)
    
    count = 0
    for i in range(0, textLen - patternLen):
        if Text[i:i+patternLen] == Pattern:
            count += 1
    
    return count


# %%    
print(PatternCount(genome[0], genome[1]))
