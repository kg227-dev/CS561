import numpy as np
from sklearn.cluster import KMeans
from spectrum_kernel import *


def parse_fasta_file(file_name):
    sequences = []
    class_names = []
    with open(file_name) as file:
        sequence = ""
        class_name = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                # Pull class name from header line
                class_name = line.split("/")[-1][6:]
                class_names.append(class_name)
                if sequence:
                    # Add sequence to list
                    sequences.append(sequence)
                    sequence = ""
            else:
                # Append sequence to current sequence string
                sequence += line
        # Add last sequence to list
        sequences.append(sequence)
    return sequences, class_names


def kmeans_clustering(sequences, kmer_size, num_clusters):
    # Generate k-mer dictionary
    kmers = set()
    for seq in sequences:
        for i in range(len(seq) - kmer_size + 1):
            kmers.add(seq[i:i+kmer_size])

    kmer_dict = {}
    for i, kmer in enumerate(sorted(kmers)):
        kmer_dict[kmer] = i

    # Convert sequences to feature vectors
    feature_vectors = np.zeros((len(sequences), len(sequences)))
    for i, seq in enumerate(sequences):
        for j, seq2 in enumerate(sequences):
            feature_vectors[i, j] = spectrum_kernel(seq, seq2, kmer_size)

    # Perform K-means clustering
    kmeans = KMeans(n_clusters=num_clusters,
                    random_state=0).fit(feature_vectors)

    # Return cluster assignments
    return kmeans.labels_


if __name__ == '__main__':
    # Parse FASTA file into sequences and class names
    sequences, class_names = parse_fasta_file(
        "Assignment2/kmeans/kmeans.fasta")

    labels = kmeans_clustering(
        sequences[0:200], 3, 5)
    print(labels)
