def make_dict(sequences, kmer_size):
    """
    Converts a list of sequences into a dictionary of kmers and their indexes.

    :param sequences: A list of DNA sequences.
    :type sequences: list of str
    :param kmer_size: The length of each k-mer.
    :type kmer_size: int
    :return: A dictionary where the keys are k-mers and the values are their indexes.
    :rtype: dict

    :author: Kush Gulati
    """
    kmers = set()
    for seq in sequences:
        for i in range(len(seq) - kmer_size + 1):
            kmers.add(seq[i:i+kmer_size])

    kmer_dict = {}
    for i, kmer in enumerate(sorted(kmers)):
        kmer_dict[kmer] = i

    return kmer_dict


def parse_fasta_file(file_name):
    """
    Converts a FASTA file to a list of sequences and a list of class names.

    :param file_name: The name of the input FASTA file.
    :type file_name: str
    :return: A tuple of two lists: the first list contains the DNA sequences in the FASTA file,
        and the second list contains the corresponding class names.
    :rtype: tuple (list of str, list of str)

    :author: Kush Gulati
    :acknowledgements: This function was written with the help of Coding4Medicine Online Textbook (https://coding4medicine.com/backup/Python/reading-fasta-files.html).
    """
    sequences = []
    class_names = []
    with open(file_name) as file:
        sequence = ""
        class_name = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                class_name = line.split("/")[-1][6:]
                class_names.append(class_name)
                if sequence:
                    sequences.append(sequence)
                    sequence = ""
            else:
                sequence += line
        sequences.append(sequence)
    return sequences, class_names