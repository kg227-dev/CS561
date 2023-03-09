import numpy as np
from spectrum_kernel import *
from extra_functions import *
import pandas as pd


def knn(X_train, y_train, X_test, k, kmer_size):
    """
    Classifies the given test sequences using the k-nearest neighbors algorithm.

    :param X_train: A list of training sequences.
    :type X_train: list
    :param y_train: A list of labels corresponding to the training sequences.
    :type y_train: list
    :param X_test: A list of test sequences to be classified.
    :type X_test: list
    :param k: The number of nearest neighbors to consider for classification.
    :type k: int
    :param kmer_size: The length of k-mers to use for feature extraction.
    :type kmer_size: int

    :return: A list of predicted labels for the test sequences.
    :rtype: list

    :author: Kush Gulati
    """
    kmer_dict = make_dict(X_train, kmer_size)
    preds = []

    X_train_feature_vectors = np.zeros((len(X_train), len(kmer_dict)))
    X_test_feature_vectors = np.zeros((len(X_test), len(kmer_dict)))

    # vectorize sequences

    for i, seq in enumerate(X_train):
        X_train_feature_vectors[i] = seq_to_vec(seq, kmer_dict, kmer_size)

    for i, seq in enumerate(X_test):
        X_test_feature_vectors[i] = seq_to_vec(seq, kmer_dict, kmer_size)

    # iterate through test vectors
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
        
    # return predicted labels
    return preds


def get_accuracy(preds, y_test):
    """
    Calculates the accuracy of a set of predictions given the correct labels.

    :param preds: A list of predicted labels.
    :type preds: list
    :param y_test: A list of correct labels.
    :type y_test: list
    :return: A float representing the proportion of correct predictions.
    :rtype: float

    :author: Kush Gulati
    """
    count = sum(1 for i in range(len(preds)) if preds[i] == y_test[i])
    accuracy = count / len(y_test)
    return accuracy


def get_output(sample_sizes, kmer_sizes, num_neighbors):
    """
    Runs the K-nearest neighbor classification algorithm with different parameter values and generates a Pandas
    DataFrame with the results.

    :param sample_sizes: A list of integers representing the sizes of the training samples to use.
    :type sample_sizes: list of int
    :param kmer_sizes: A list of integers representing the lengths of the k-mers to use.
    :type kmer_sizes: list of int
    :param num_neighbors: A list of integers representing the number of neighbors to consider.
    :type num_neighbors: list of int
    :return: A Pandas DataFrame with the accuracy of the K-nearest neighbor algorithm for each combination of parameters.
    :rtype: pandas.DataFrame

    :author: Kush Gulati
    """
    data = {'sample_size': [],
            'kmer_size': [],
            'num_neighbors': [],
            'accuracy': []}
    df = pd.DataFrame(data)
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
                accuracy = get_accuracy(preds, y_test)
                new_row = {'sample_size': sample_size, 'kmer_size': kmer_size,
                           'num_neighbors': num_neighbor, 'accuracy': accuracy}
                df = df.append(new_row, ignore_index=True)
                print("Calculated accuracy for sample size {}, kmer size {}, and # neighbors {}.".format(
                    sample_size, kmer_size, num_neighbor))
    return df


if __name__ == "__main__":
    sample_sizes = ["10", "30", "100"]
    kmer_sizes = [2, 4, 6, 8]
    num_neighbors = [1, 3, 5, 7]
    df = get_output(sample_sizes, kmer_sizes, num_neighbors)
    df.to_csv('knn_metrics.csv', index=False)
