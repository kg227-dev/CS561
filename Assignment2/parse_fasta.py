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