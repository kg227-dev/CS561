import numpy as np
from parse_fasta import *
from spectrum_kernel import *
import pandas as pd


def make_dict(sequences, kmer_size):
    kmers = set()
    for seq in sequences:
        for i in range(len(seq) - kmer_size + 1):
            kmers.add(seq[i:i+kmer_size])

    kmer_dict = {}
    for i, kmer in enumerate(sorted(kmers)):
        kmer_dict[kmer] = i

    return kmer_dict


def knn(X_train, y_train, X_test, k, kmer_size):
    kmer_dict = make_dict(X_train, kmer_size)
    preds = []

    X_train_feature_vectors = np.zeros((len(X_train), len(kmer_dict)))
    X_test_feature_vectors = np.zeros((len(X_test), len(kmer_dict)))

    for i, seq in enumerate(X_train):
        X_train_feature_vectors[i] = seq_to_vec(seq, kmer_dict, kmer_size)

    for i, seq in enumerate(X_test):
        X_test_feature_vectors[i] = seq_to_vec(seq, kmer_dict, kmer_size)

    for i in range(len(X_test_feature_vectors)):
        similarities = []
        for j in range(len(X_train_feature_vectors)):
            similarity = spectrum_kernel(
                X_test_feature_vectors[i], X_train_feature_vectors[j], kmer_size, kmer_dict)
            similarities.append((similarity, y_train[j]))

        sorted_similarities = sorted(similarities, reverse=True)
        k_nearest_neighbors = sorted_similarities[:k]

        counts = {}
        for neighbor in k_nearest_neighbors:
            label = neighbor[1]
            if label in counts:
                counts[label] += 1
            else:
                counts[label] = 1

        max_count = 0
        max_label = None
        for label, count in counts.items():
            if count > max_count:
                max_count = count
                max_label = label
        preds.append(max_label)
    return preds


def get_accuracy(sample_size, kmer_size, num_neighbor, preds, y_test):
    # print("\nSample Size: {}".format(sample_size))
    # print("Kmer Size: {}".format(kmer_size))
    # print("Number of Neighbors: {}".format(num_neighbor))
    count = sum(1 for i in range(len(preds)) if preds[i] == y_test[i])
    accuracy = count / len(y_test)
    return accuracy


if __name__ == "__main__":
    sample_size = "10"
    X_train_exon, y_train_exon = parse_fasta_file(
        "Assignment2/KNN/train-exons"+sample_size+".fasta")
    X_train_intron, y_train_intron = parse_fasta_file(
        "Assignment2/KNN/train-introns"+sample_size+".fasta")
    X_train = X_train_exon + X_train_intron
    y_train = y_train_exon + y_train_intron
    X_test, y_test = parse_fasta_file("Assignment2/KNN/test.fasta")

    data = {'sample_size': [],
            'kmer_size': [],
            'num_neighbors': [],
            'accuracy': []}
    df = pd.DataFrame(data)

    sample_sizes = ["10", "30", "100"]
    kmer_sizes = [2, 4, 6, 8]
    num_neighbors = [1, 3, 5, 7]
    for sample_size in sample_sizes:
        X_train_exon, y_train_exon = parse_fasta_file(
            "Assignment2/KNN/train-exons"+sample_size+".fasta")
        X_train_intron, y_train_intron = parse_fasta_file(
            "Assignment2/KNN/train-introns"+sample_size+".fasta")
        X_train = X_train_exon + X_train_intron
        y_train = y_train_exon + y_train_intron
        X_test, y_test = parse_fasta_file("Assignment2/KNN/test.fasta")
        for kmer_size in kmer_sizes:
            for num_neighbor in num_neighbors:
                preds = knn(X_train, y_train, X_test, num_neighbor, kmer_size)
                accuracy = get_accuracy(
                    sample_size, kmer_size, num_neighbor, preds, y_test)
                new_row = {'sample_size': sample_size, 'kmer_size': kmer_size,
                           'num_neighbors': num_neighbor, 'accuracy': accuracy}
                df = df.append(new_row, ignore_index=True)

    df.to_csv('knn_metrics.csv', index=False)
