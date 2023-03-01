import numpy as np
from spectrum_kernel import *
from sklearn.cluster import KMeans


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


def compute_similarity_matrix(sequences, k):
    n = len(sequences)
    sim_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            sim_matrix[i][j] = sim_matrix[j][i] = spectrum_kernel(sequences[i], sequences[j], k)
        print(i)
    return sim_matrix



def cluster_sequences(sequences, k, n_clusters):
    sim_matrix = compute_similarity_matrix(sequences, k)
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(sim_matrix)
    labels = kmeans.labels_
    return labels


if __name__ == '__main__':
    # Parse FASTA file into sequences and class names
    sequences, class_names = parse_fasta_file(
        "Assignment2/kmeans/kmeans.fasta")


    labels = cluster_sequences(
        sequences, 3, 4)
    print(labels)
