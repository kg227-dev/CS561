from spectrum_kernel import *
from mismatch_kernel import *


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


def make_dict(sequences, kmer_size):
    kmers = set()
    for seq in sequences:
        for i in range(len(seq) - kmer_size + 1):
            kmers.add(seq[i:i+kmer_size])

    kmer_dict = {}
    for i, kmer in enumerate(sorted(kmers)):
        kmer_dict[kmer] = i

    return kmer_dict

def kmeans_mismatch_cluster(sequences, kmer_size, num_clusters, max_iter=1000):
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
                similarities[j] = np.dot(seq, centroid)
            cluster_assignments[i] = np.argmax(similarities)

        # Update centroids
        for j in range(num_clusters):
            cluster_indices = np.where(cluster_assignments == j)[0]
            if len(cluster_indices) > 0:
                centroids[j] = np.mean(
                    feature_vectors[cluster_indices, :], axis=0)

    return cluster_assignments, centroids

def kmeans_spectrum_cluster(sequences, kmer_size, num_clusters, max_iter=100):
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
                similarities[j] = np.dot(seq, centroid)
            cluster_assignments[i] = np.argmax(similarities)

        # Update centroids
        for j in range(num_clusters):
            cluster_indices = np.where(cluster_assignments == j)[0]
            if len(cluster_indices) > 0:
                centroids[j] = np.mean(
                    feature_vectors[cluster_indices, :], axis=0)

    return cluster_assignments, centroids


def get_output(cluster_assignments, kmer_size, num_clusters):
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
    kmer_size = 3
    num_clusters = 3
    cluster_assignments, centroids = kmeans_spectrum_cluster(
        sequences, kmer_size, num_clusters)
    get_output(cluster_assignments, kmer_size, num_clusters)
