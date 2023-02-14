from parsing import *
import numpy as np

"""
Affine scoring to penalize small gaps over large gaps
By Sydney Ballard
"""
def affine(sequence1, sequence2, matrix):
    gap_open = -4
    gap_extend = -0.1
    m, n = len(sequence1), len(sequence2)

    # Get dictionary of nucletide match/mismatch score
    scoring_matrix = get_scoring_matrix_dict(matrix)

    # Initialize 3 matrices for traceback
    M = np.zeros((m+1, n+1)) # This is the main "scorekeeper", tracks only match/mismatch
    seq1_GAP = np.zeros((m+1, n+1)) # Tracks gap creation/elongation in sequence 1
    seq2_GAP = np.zeros((m+1, n+1)) # Tracks gap creation/elongation in sequence 2

    # Fill respective traceback matrices
    for i in range(1, m + 1):
        M[i][0] = gap_open + i * gap_extend
        seq1_GAP[i][0] = M[i][0] + gap_extend
        seq2_GAP[i][0] = float("-inf")
        
    for j in range(1, n + 1):
        M[0][j] = gap_open + j * gap_extend
        seq1_GAP[0][j] = float("-inf")
        seq2_GAP[0][j] = M[0][j] + gap_extend
        
    # Calculate matrix cell scores based on scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Update score from scoring matrix
            match_score = M[i - 1][j - 1] + scoring_matrix[sequence1[i - 1] + sequence2[j - 1]] 

            ### Fill out next matrix cell for each matrix ###

            # Matrix tracking best alignment of gap (creating one or extending one) in sequence 1
            seq1_GAP[i][j] = max(seq1_GAP[i - 1][j] + gap_extend, M[i - 1][j] + gap_open + gap_extend)
            
            # Matrix tracking best alignment of gap (creating one or extending one) in sequence 2
            seq2_GAP[i][j] = max(seq2_GAP[i][j - 1] + gap_extend, M[i][j - 1] + gap_open + gap_extend)

            # Matrix tracking best alignment--will add/extend gap or keep match/mismatch based on previous calculations
            M[i][j] = max(match_score, seq1_GAP[i][j], seq2_GAP[i][j])

    # Optimal alignment can be found from backtracking from M[m][n]
    align1, align2 = [], []
    i, j = m, n
    
    # Exhaustive backtrace
    while i > 0 and j > 0:
        if M[i][j] == seq1_GAP[i][j]:
            align1.append(sequence1[i - 1])
            align2.append("-")
            i -= 1
        elif M[i][j] == seq2_GAP[i][j]:
            align1.append("-")
            align2.append(sequence2[j - 1])
            j -= 1
        else:
            align1.append(sequence1[i - 1])
            align2.append(sequence2[j - 1])
            i -= 1
            j -= 1
        
    while i > 0:
        align1.append(sequence1[i - 1])
        align2.append("-")
        i -= 1
        
    while j > 0:
        align1.append("-")
        align2.append(sequence2[j - 1])
        j -= 1

    sequence_alignment = "".join(align1[::-1]), "".join(align2[::-1])
        
    return sequence_alignment