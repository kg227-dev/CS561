import random


def mismatch_kernel(kmers1, kmers2, k):
    # Compute all k-mers for both sequences
    kmers1 = [seq1[i:i+k] for i in range(len(seq1)-k+1)]
    kmers2 = [seq2[i:i+k] for i in range(len(seq2)-k+1)]

    # Initialize similarity score to 0
    similarity = 0


    # Loop over all k-mers in first sequence
    for kmer1 in kmers1:
        # Loop over all k-mers in second sequence
        for kmer2 in kmers2:
            # Count number of mismatches
            mismatches = sum([1 for i in range(k) if kmer1[i] != kmer2[i]])

            # Update the similarity score based on the number of mismatches
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
