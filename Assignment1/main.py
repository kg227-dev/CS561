import pandas as pd
from parsing import *
from affine import *
from linear import *
from get_alignment_output import *


def main():
    ############################
    # PAIRWISE SCORING MATRIX  #
    ############################
    # File path to scoring matrix
    matrix_file = "Assignment1/input files/matrix.txt"

    ############################
    #   SELECTING FILE PATHS   #
    ############################
    # File paths ran two at a time:
    #   - close-first.fasta  + close-second.fasta
    #   - distant-first.fast + distant-second.fasta

    # CLOSE
    first_file = "Assignment1/input files/close-first.fasta"
    second_file = "Assignment1/input files/close-second.fasta"

    # DISTANT
    # first_file = "Assignment1/input files/distant-first.fasta"
    # second_file = "Assignment1/input files/distant-second.fasta"

    ############################
    # SCORING / ALIGNMENT TYPE #
    ############################
    # output() function performs bulk of operations, including console output to new file

    # Linear scoring:
    #output(first_file, second_file, "linear-close.txt", matrix_file)

    # Affine scoring:
    output(first_file, second_file, "affine-distant.txt", matrix_file, 1)


if __name__ == "__main__":
    main()
