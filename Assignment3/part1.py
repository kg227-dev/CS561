

def one_hot_encoding(dna_sequence):
    """
    Converts a DNA sequence to one-hot-encoding

    :param dna_sequence: nucleotide sequence
    :type dna_sequence: str
    :return: lst of encoding for each nucleotide in the original DNA string
    :rtype: lst of lst of ints

    :author: Sydney Ballard
    """
    # Define a dictionary for 1-hot encoding
    nucleotide_dict = {'A': [1, 0, 0, 0], 
                       'C': [0, 1, 0, 0], 
                       'G': [0, 0, 1, 0], 
                       'T': [0, 0, 0, 1]}
    
    # Convert the DNA sequence to 1-hot encoding
    one_hot_encoding = [nucleotide_dict[base] for base in dna_sequence]  

    # Return list of one hot encodings
    return one_hot_encoding


def parse_data(file_name):
    """
    Converts a FASTA file to a dictionary.

    :param file_name: The name of the input FASTA file.
    :type file_name: str
    :return: A dictionary of each key with a list of two elements: 
        the key is the original nucleotide sequence
        the first element of the list contains the one-hot-encoding of the sequence (reference one_hot_encoding function above)
        and the second element of the list contains the corresponding class names/enhancer identifier (0 or 1)
    :rtype: dict

    :author: Sydney Ballard
    :acknowledgements: Adapted from parse_fasta_file, which Kush Gulati wrote for previous assignments
    """
    # Initialize dictionary 
    # {sequence: [one hot encoding of sequence, enhancer identifier]}
    sequence_data = {} # dict will hold 0/1 for class and sequence

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
    
    for idx in range(len(sequences)):
        seq = sequences[idx]
        class_encoding = class_names[idx]

        # 0th element = encoding, 1st element = enhancer identifier
        sequence_data[seq] = [one_hot_encoding(seq), class_encoding]

    # print(list(sequence_data.values())[0][1])    
    return sequence_data