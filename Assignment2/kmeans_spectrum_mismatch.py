from spectrum_kernel import *
from mismatch_kernel import *
from extra_functions import *


def kmeans_cluster(sequences, kernel_type, kmer_size, num_clusters, max_iter=100):
    """
    Performs k-means clustering on a list of sequences using the appropriate kernel.

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

    :author: Kush Gulati & Sydney Ballard
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

                ### appropriate kernel operation chosen
                if kernel_type == "SPECTRUM":
                    kernel_result = spectrum_kernel(
                        seq, centroid, kmer_size, kmer_dict)
                elif kernel_type == "MISMATCH":
                    kernel_result = mismatch_kernel(
                        seq, centroid, kmer_size)
                    
                # assign appropriate score to matrix
                similarities[j] = kernel_result
                # if kernel_result != 0: print(kernel_result) 
                    
            cluster_assignments[i] = np.argmax(similarities)


        for j in range(num_clusters):
            cluster_indices = np.where(cluster_assignments == j)[0]
            if len(cluster_indices) > 0:
                centroids[j] = np.mean(
                    feature_vectors[cluster_indices, :], axis=0)

    return cluster_assignments, centroids


def get_output(kernel_type, cluster_assignments, kmer_size, num_clusters):
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

    print("{} KERNEL (KMER={}, CLUSTERS={}):".format(kernel_type, kmer_size, num_clusters))
    for i in range(num_clusters):
        print("CLUSTER {}:".format(i+1))
        indices = np.where(cluster_assignments == i)[0]
        classes = [class_names[index] for index in indices]
        freq_dict = {s: classes.count(s)
                     for s in ['intergenic', 'intron', 'exon']}
        
        intergenic_count = freq_dict['intergenic']
        intron_count = freq_dict['intron']
        exon_count = freq_dict['exon']

        print("\tintergenic = {} ({})".format(avg(intergenic_count, len(classes)), intergenic_count))
        print("\tintron = {} ({})".format(avg(intron_count, len(classes)), intron_count))
        print("\texon = {} ({})".format(avg(exon_count, len(classes)), exon_count))
    return None

# implemented function to find average, rounded to next hundredth
# bypasses divide by zero issue when using mismatch kernel
def avg(num, divide_by):
    if divide_by != 0:
        return np.round(num/divide_by, 2)
    else:
        return num


if __name__ == '__main__':
    # Parse FASTA file into sequences and class names
    sequences, class_names = parse_fasta_file(
        "/Users/sydneyballard/Desktop/Desktop - Sydney’s MacBook Pro/CS 561/cs561 repository COLLABORATIVE W KUSH/CS561/Assignment2/kmeans/kmeans.fasta")
    # kernel = "MISMATCH"
    kernel = "SPECTRUM"
    kmer_size = 2
    num_clusters = 2
    cluster_assignments, centroids = kmeans_cluster(
        sequences, kernel, kmer_size, num_clusters)
    
    print(cluster_assignments)
    print(centroids)

    get_output(kernel, cluster_assignments, kmer_size, num_clusters)
