from spectrum_kernel import *
from extra_functions import *


def kmeans_cluster(sequences, kmer_size, num_clusters, max_iter=100):
    """
    Performs k-means clustering on a list of sequences using the spectrum kernel.

    :param sequences: A list of input sequences.
    :type sequences: list
    :param kmer_size: The length of k-mers to use for feature extraction.
    :type kmer_size: int
    :param num_clusters: The number of clusters to form.
    :type num_clusters: int
    :param max_iter: The maximum number of iterations to perform.
    :type max_iter: int

    :return: A tuple containing a list of cluster assignments for each input sequence and a list of centroid sequences.
    :rtype: tuple

    :author: Kush Gulati
    """
    # Generate k-mer dictionary
    kmer_dict = make_dict(sequences, kmer_size)

    # Initialize centroids
    centroids = np.random.choice(sequences, size=num_clusters, replace=False)

    centroids = [seq_to_vec(centroid, kmer_dict, kmer_size)
                 for centroid in centroids]

    # Convert sequences to feature vectors
    feature_vectors = np.zeros((len(sequences), len(kmer_dict)))
    for i, seq in enumerate(sequences):
        feature_vectors[i, :] = seq_to_vec(seq, kmer_dict, kmer_size)

    cluster_assignments = np.zeros(len(sequences))

    # K-means algorithm
    for iter in range(max_iter):
        # Assign sequences to closest centroid
        for i, seq in enumerate(feature_vectors):
            similarities = np.zeros(num_clusters)
            for j, centroid in enumerate(centroids):
                similarities[j] = spectrum_kernel(
                    seq, centroid, kmer_size, kmer_dict)
            cluster_assignments[i] = np.argmax(similarities)

        # Update centroids
        for j in range(num_clusters):
            cluster_indices = np.where(cluster_assignments == j)[0]
            if len(cluster_indices) > 0:
                centroids[j] = np.mean(
                    feature_vectors[cluster_indices, :], axis=0)

    return cluster_assignments, centroids


def get_output(cluster_assignments, kmer_size, num_clusters):
    """
    Print the output of clustering analysis with spectrum kernel.

    :param cluster_assignments: numpy array containing the cluster assignments 
    for each sequence
    :type cluster_assignments: numpy.ndarray
    :param kmer_size: size of k-mers used in the spectrum kernel
    :type kmer_size: int
    :param num_clusters: number of clusters used in the K-means algorithm
    :type num_clusters: int
    :return: None

    :author: Kush Gulati
    """

    print("SPECTRUM KERNEL (KMER={}, CLUSTERS={}):".format(kmer_size, num_clusters))
    for i in range(num_clusters):
        print("CLUSTER {}:".format(i+1))
        indices = np.where(cluster_assignments == i)[0]
        classes = [class_names[index] for index in indices]
        freq_dict = {s: classes.count(s)
                     for s in ['intergenic', 'intron', 'exon']}
        print("\tintergenic = {} ({})".format(
            np.round(freq_dict['intergenic']/len(classes), 2), freq_dict['intergenic']))
        print("\tintron = {} ({})".format(
            np.round(freq_dict['intron']/len(classes), 2), freq_dict['intron']))
        print("\texon = {} ({})".format(
            np.round(freq_dict['exon']/len(classes), 2), freq_dict['exon']))
    return None


if __name__ == '__main__':
    # Parse FASTA file into sequences and class names
    sequences, class_names = parse_fasta_file(
        "Assignment2/kmeans/kmeans.fasta")
    kmer_size = 6
    num_clusters = 5
    cluster_assignments, centroids = kmeans_cluster(
        sequences, kmer_size, num_clusters)
    get_output(cluster_assignments, kmer_size, num_clusters)
