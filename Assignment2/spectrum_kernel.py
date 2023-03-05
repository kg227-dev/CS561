import numpy as np
import random


def spectrum_kernel(seq1, seq2, k, kmer_dict=None):
    if kmer_dict is None:
        # Convert sequences to feature vectors
        vec1 = seq_to_vec(seq1, kmer_dict, k)
        vec2 = seq_to_vec(seq2, kmer_dict, k)
    # Compute dot product similarity
        return np.dot(vec1, vec2)
    else:
        return np.dot(seq1, seq2)


def seq_to_vec(seq, kmer_dict, k):
    # Initialize feature vector
    vec = np.zeros(len(kmer_dict))

    # Populate feature vector with k-mer counts
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in kmer_dict:
            vec[kmer_dict[kmer]] += 1

    return vec


if __name__ == "__main__":
    nucleotides = ['A', 'C', 'G', 'T']
    seq1 = ''.join(random.choice(nucleotides) for i in range(100))
    seq2 = ''.join(random.choice(nucleotides) for i in range(100))

    k = 3
    similarity = spectrum_kernel(seq1, seq2, k)

    print("Sequence 1:", seq1)
    print("Sequence 2:", seq2)
    print("K-mer size:", k)
    print("Similarity:", similarity)
