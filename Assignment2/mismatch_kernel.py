import random


"""
# the mismatch kernel is more appropriate for comparing sequences 
# that may have mutations or errors
"""


def mismatch_kernel(seq1, seq2, k):
    """
    Compute the mismatch kernel similarity_score score between two DNA sequences based on k-mers.

    :param seq1: A DNA sequence represented as a string.
    :param seq2: A DNA sequence represented as a string.
    :param k: The length of the k-mers to be considered.
    :return: An integer representing the similarity_score score between the two sequences based on the mismatch kernel.

    :author: Kush Gulati + Sydney Ballard
    """
    n = len(seq1)
    kernel = 0
    
    for i in range(n-k+1):
        kmer1 = seq1[i:i+k]
        kmer2 = seq2[i:i+k]
        
        if hamming_distance(kmer1, kmer2) <= 1:
            kernel += 1
            
    return float(kernel)


def hamming_distance(str1, str2):
    """
    Computes the Hamming distance between two strings of equal length.
    
    :param str1: the first string
    :param str2: the second string
    :return: the Hamming distance between the two strings
    """
    assert len(str1) == len(str2)
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))


if __name__ == "__main__":
    nucleotides = ['A', 'C', 'G', 'T']

    # randomly creating an nucleotide sequence for testing
    seq1 = ''.join(random.choice(nucleotides) for i in range(100))
    seq2 = ''.join(random.choice(nucleotides) for i in range(100))

    k = 3
    similarity_score = mismatch_kernel(seq1, seq2, k)

    print("Sequence 1: ", seq1)
    print("Sequence 2: ", seq2)
    print("similarity_score score (k={}): {}".format(k, similarity_score))
