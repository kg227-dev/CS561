import numpy as np
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
    # Compute pairwise similarities between sequences
    n = len(sequences)
    similarities = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            sim = spectrum_kernel(sequences[i], sequences[j], kmer_size)
            similarities[i, j] = sim
            similarities[j, i] = sim

    # Initialize cluster centroids randomly
    centroids = np.random.rand(num_clusters, n)

    # Iterate until convergence
    for i in range(10):
        # Assign each sequence to nearest centroid
        distances = np.linalg.norm(similarities[:, :, np.newaxis] - centroids.T[np.newaxis, np.newaxis, :, :], axis=3)
        cluster_indices = np.argmin(distances, axis=2)

        # Update centroids
        feature_vectors = np.zeros((num_clusters, n))
        for cluster in range(num_clusters):
            feature_vectors[cluster, :] = np.sum(similarities[cluster_indices == cluster, :], axis=0)
        nonempty_clusters = np.sum(cluster_indices == np.arange(num_clusters)[:, np.newaxis], axis=1)
        centroids[nonempty_clusters > 0, :] = feature_vectors[nonempty_clusters > 0, :] / nonempty_clusters[nonempty_clusters > 0, np.newaxis]

    # Return cluster assignments
    return np.argmax(distances, axis=1)

if __name__ == '__main__':
    # Parse FASTA file into sequences and class names
    sequences, class_names = parse_fasta_file(
        "Assignment2/kmeans/kmeans.fasta")


    labels = kmeans_clustering(
        sequences[0:10], 3, 4)
