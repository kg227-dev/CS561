import numpy as np
import random


def spectrum_kernel(seq1, seq2, k, kmer_dict=None):
    """
    Compute the spectrum kernel similarity score between two DNA sequences.

    :param seq1: The first DNA sequence, either as a string or as a feature vector.
    :type seq1: Union[str, np.ndarray]
    :param seq2: The second DNA sequence, either as a string or as a feature vector.
    :type seq2: Union[str, np.ndarray]
    :param k: The length of the k-mers.
    :type k: int
    :param kmer_dict: A dictionary mapping k-mers to indices, used for vectorization.
    :type kmer_dict: Optional[Dict[str, int]]
    :return: The similarity score between the two sequences.
    :rtype: float

    :author: Kush Gulati
    """
    if kmer_dict is None:
        kmers = set()
        for i in range(len(seq1) - k + 1):
            kmers.add(seq1[i:i+k])
        for i in range(len(seq2) - k + 1):
            kmers.add(seq2[i:i+k])
        kmer_dict = {}
        for i, kmer in enumerate(sorted(kmers)):
            kmer_dict[kmer] = i
        vec1 = seq_to_vec(seq1, kmer_dict, k)
        vec2 = seq_to_vec(seq2, kmer_dict, k)
        return np.dot(vec1, vec2)
    else:
        return np.dot(seq1, seq2)


def seq_to_vec(seq, kmer_dict, k):
    """
    Converts a DNA sequence into a k-mer frequency vector using a pre-defined k-mer dictionary.

    :param seq: A string representing a DNA sequence.
    :type seq: str
    :param kmer_dict: A dictionary where keys are k-mers and values are their corresponding indices in the output vector.
    :type kmer_dict: dict
    :param k: An integer representing the k-mer size.
    :type k: int
    :return: A numpy array representing the k-mer frequency vector of the input sequence.
    :rtype: np.ndarray

    :author: Kush Gulati
    """
    vec = np.zeros(len(kmer_dict))

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
