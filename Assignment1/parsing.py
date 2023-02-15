import re
import pandas as pd

# function that reads FASTA file and returns a dictionary
# https://coding4medicine.com/backup/Python/reading-fasta-files.html


def fasta(file):
    f = open(file, 'r')
    lines = f.readlines()
    hre = re.compile('>(\S+)')
    lre = re.compile('^(\S+)$')
    gene = {}
    for line in lines:
        outh = hre.search(line)
        if outh:
            id = outh.group(1)
        else:
            outl = lre.search(line)
            if (id in gene.keys()):
                gene[id] += outl.group(1)
            else:
                gene[id] = outl.group(1)
    return gene


"""
Function returns the matrix.txt file as a dictionary in the form { pairing : score }
By: Sydney Ballard
"""


def get_scoring_matrix_dict(matrix_txt):
    matrix_df = pd.read_csv(matrix_txt, sep=' ',
                            skipinitialspace=True, index_col=0)
    extracted_matrix = [[x, y, matrix_df[x][y]]
                        for x in matrix_df.index for y in matrix_df.columns]

    scoring_dict = {}

    for score in extracted_matrix:
        key = score[0] + score[1]
        if key not in scoring_dict:
            scoring_dict[key] = score[2]

    return scoring_dict
