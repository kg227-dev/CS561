import random


def mismatch_kernel(seq1, seq2, k):
    """
    Compute the mismatch kernel similarity score between two DNA sequences based on k-mers.

    :param seq1: A DNA sequence represented as a string.
    :param seq2: A DNA sequence represented as a string.
    :param k: The length of the k-mers to be considered.
    :return: An integer representing the similarity score between the two sequences based on the mismatch kernel.

    :author: Kush Gulati
    """

    kmers1 = [seq1[i:i+k] for i in range(len(seq1)-k+1)]
    kmers2 = [seq2[i:i+k] for i in range(len(seq2)-k+1)]

    similarity = 0

    for kmer1 in kmers1:
        for kmer2 in kmers2:
            mismatches = sum([1 for i in range(k) if kmer1[i] != kmer2[i]])

            if mismatches <= 1:
                similarity += 1

    return similarity


if __name__ == "__main__":
    nucleotides = ['A', 'C', 'G', 'T']
    seq1 = ''.join(random.choice(nucleotides) for i in range(5))
    seq2 = ''.join(random.choice(nucleotides) for i in range(5))

    k = 3
    similarity = mismatch_kernel(seq1, seq2, k)

    print("Sequence 1: ", seq1)
    print("Sequence 2: ", seq2)
    print("Similarity score (k={}): {}".format(k, similarity))
