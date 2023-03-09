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

    :author: Kush Gulati
    """

    seq1_kmers = [seq1[i:i+k] for i in range(len(seq1)-k+1)] # sequence 1 list of kmers of length k
    seq2_kmers = [seq2[i:i+k] for i in range(len(seq2)-k+1)] # sequence 2 list of kmers of length k

    similarity_score = 0 # initializing similarity_score score to 0

    ## The nested loops below compare each respective kmer

    # loop through each kmer of seq1_kmers list
    for kmer1 in seq1_kmers: 
        # loop through each kmer of seq2_kmers list
        for kmer2 in seq2_kmers:
            # increment mismatch score by 1 if mismatch if found
            mismatches = sum([1 for i in range(k) if kmer1[i] != kmer2[i]])

            # if mismatches <= 1:
            #     similarity_score += 1

            similarity_score += 2**(-mismatches)

    return similarity_score


if __name__ == "__main__":
    nucleotides = ['A', 'C', 'G', 'T']

    # randomly creating an nucleotide sequence for testing
    seq1 = ''.join(random.choice(nucleotides) for i in range(5))
    seq2 = ''.join(random.choice(nucleotides) for i in range(5))

    k = 3
    similarity_score = mismatch_kernel(seq1, seq2, k)

    print("Sequence 1: ", seq1)
    print("Sequence 2: ", seq2)
    print("similarity_score score (k={}): {}".format(k, similarity_score))
