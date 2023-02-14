import pandas as pd
from parsing import *
from affine import *
from linear import *
from get_alignment_output import *

def main():
    matrix_file ="cs561 repository COLLABORATIVE W KUSH/CS561/Assignment1/matrix.txt"


    ### CLOSE
    # first_file = "cs561 repository COLLABORATIVE W KUSH/CS561/Assignment1/close-first.fasta"
    # second_file = "cs561 repository COLLABORATIVE W KUSH/CS561/Assignment1/close-second.fasta"

    ### DISTANT
    first_file = "cs561 repository COLLABORATIVE W KUSH/CS561/Assignment1/distant-first.fasta"
    second_file = "cs561 repository COLLABORATIVE W KUSH/CS561/Assignment1/distant-second.fasta"

    # Linear scoring: 
    # output(first_file, second_file, "linear-close.txt", matrix_file)

    # Affine scoring:
    output(first_file, second_file, "affine-distant.txt", matrix_file, 1)
    


if __name__ == "__main__":
    main()