import numpy as np

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

# %%
data = np.loadtxt('./data/dataset_9_3.txt', dtype='str')
print('p: {0}\n\nq: {1}\n\nHamming Distance: {2}'.format(data[0], data[1], HammingDistance(data[0], data[1])))