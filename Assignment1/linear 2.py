import numpy as np
import math
from parsing import *

"""
Function returns the score of each pairing based on matrix.txt
By: Kush Gulati
"""


def get_score(a, b):
    if a == b:
        score = 5
    elif ((a == 'A') and (b == 'G')) or ((a == 'G') and (b == 'A')):
        score = -1
    elif ((a == 'C') and (b == 'T')) or ((a == 'T') and (b == 'C')):
        score = -1
    elif ((a == 'A') and (b == 'T')) or ((a == 'T') and (b == 'A')):
        score = -4
    elif ((a == 'C') and (b == 'G')) or ((a == 'G') and (b == 'C')):
        score = -4
    elif a == '-' and b == '-':
        score = -1
    else:
        score = -2
    return score


"""
Implementation of Needleman_Wunsch algorithm with linear scoring
By: Kush Gulati
"""


def needleman_wunsch(x, y):
    # find length of both sequences
    lenx = len(x)
    leny = len(y)

    # create matrix of zeros
    F = np.zeros((leny+1, lenx+1))

    # set gap penalty
    gap_pen = -1

    # fill first column
    for i in range(0, leny+1):
        F[i][0] = gap_pen * i

    # fill first row
    for j in range(0, lenx+1):
        F[0][j] = gap_pen * j

    # fill values in score matrix
    for i in range(1, len(y)+1):
        for j in range(1, len(x)+1):
            # check top, left, and diag cells
            match = F[i-1, j-1] + get_score(x[j-1], y[i-1])
            deletion = F[i-1][j] + gap_pen
            insertion = F[i][j-1] + gap_pen
            # get the max score
            F[i, j] = max(match, deletion, insertion)

    score = F[i, j]

    # store the alignment
    linex = ""
    liney = ""

    # set position of i and j to bottom right of matrix
    i = leny
    j = lenx

    # end traceback when matrix reaches the top or the left side
    while i > 0 and j > 0:
        current = F[i][j]
        diag = F[i-1][j-1]
        up = F[i][j-1]
        left = F[i-1][j]
        if current == diag + get_score(x[j-1], y[i-1]):
            linex += x[j-1]
            liney += y[i-1]
            i -= 1
            j -= 1
        elif current == up + gap_pen:
            linex += x[j-1]
            liney += '-'
            j -= 1
        elif current == left + gap_pen:
            linex += '-'
            liney += y[i-1]
            i -= 1

    # trace up to the top left cell
    while j > 0:
        linex += x[j-1]
        liney += '-'
        j -= 1
    while i > 0:
        linex += '-'
        liney += y[i-1]
        i -= 1

    # traverse the score matrix from the bottom right, so reverse the two sequences
    linex = linex[::-1]
    liney = liney[::-1]

    return (linex, liney, score)


"""
Function that returns the mean indel length of a given alignment
By: Kush Gulati
"""


def get_mean_indel_length(s):
    lengths = []
    current_length = 0
    for char in s:
        if char == " ":
            current_length += 1
        elif current_length > 0:
            lengths.append(current_length)
            current_length = 0
    if current_length > 0:
        lengths.append(current_length)
    return sum(lengths)/len(lengths)


"""
Function that properly formats the output and prints it to the console
By: Kush Gulati
"""


def get_output(x, y):
    linex, liney, score = needleman_wunsch(x, y)
    out_line = ""
    for i in range(max(len(linex), len(liney))):
        if linex[i] == liney[i]:
            out_line += "|"
        elif (linex[i] == "-" or liney[i] == "-"):
            out_line += " "
        else:
            out_line += "*"
    matches = out_line.count('|')
    indels = out_line.count(' ')
    percent_identity = round(100 * matches / ((len(linex)+len(liney))/2))
    print("Matches: " + str(matches))
    print("Percent identity: " + str(percent_identity)+"%")
    print("Indels: " + "number="+str(indels) + ", mean length=" +
          str(round(get_mean_indel_length(out_line), 1)))
    print("Alignment length: " + str(max(len(linex), len(liney))))
    print("Score="+str(score))
    for i in range(math.ceil(max(len(linex), len(liney))/60)):
        print("\n"+linex[60*i:60*(i+1)])
        print(out_line[60*i:60*(i+1)])
        print(liney[60*i:60*(i+1)])


"""
Linear-close
x = fasta("Assignment1/input files/close-first.fasta").get("first1")
y = fasta("Assignment1/input files/close-second.fasta").get("second1")
print("\nAlignment 1:")
print("\nSequence 1: close-first")
print("Sequence 2: close-second")
get_output(x,y)

x = fasta("Assignment1/input files/close-first.fasta").get("first2")
y = fasta("Assignment1/input files/close-second.fasta").get("second2")
print("\nAlignment 2:")
print("\nSequence 1: close-first")
print("Sequence 2: close-second")
get_output(x,y)

x = fasta("Assignment1/input files/close-first.fasta").get("first3")
y = fasta("Assignment1/input files/close-second.fasta").get("second3")
print("\nAlignment 3:")
print("\nSequence 1: close-first")
print("Sequence 2: close-second")
get_output(x,y)

x = fasta("Assignment1/input files/close-first.fasta").get("first4")
y = fasta("Assignment1/input files/close-second.fasta").get("second4")
print("\nAlignment 4:")
print("\nSequence 1: close-first")
print("Sequence 2: close-second")
get_output(x,y)

x = fasta("Assignment1/input files/close-first.fasta").get("first5")
y = fasta("Assignment1/input files/close-second.fasta").get("second5")
print("\nAlignment 5:")
print("\nSequence 1: close-first")
print("Sequence 2: close-second")
get_output(x,y)

x = fasta("Assignment1/input files/close-first.fasta").get("first6")
y = fasta("Assignment1/input files/close-second.fasta").get("second6")
print("\nAlignment 6:")
print("\nSequence 1: close-first")
print("Sequence 2: close-second")
get_output(x,y)

x = fasta("Assignment1/input files/close-first.fasta").get("first7")
y = fasta("Assignment1/input files/close-second.fasta").get("second7")
print("\nAlignment 7:")
print("\nSequence 1: close-first")
print("Sequence 2: close-second")
get_output(x,y)

x = fasta("Assignment1/input files/close-first.fasta").get("first8")
y = fasta("Assignment1/input files/close-second.fasta").get("second8")
print("\nAlignment 8:")
print("\nSequence 1: close-first")
print("Sequence 2: close-second")
get_output(x,y)

x = fasta("Assignment1/input files/close-first.fasta").get("first9")
y = fasta("Assignment1/input files/close-second.fasta").get("second9")
print("\nAlignment 9:")
print("\nSequence 1: close-first")
print("Sequence 2: close-second")
get_output(x,y)

x = fasta("Assignment1/input files/close-first.fasta").get("first10")
y = fasta("Assignment1/input files/close-second.fasta").get("second10")
print("\nAlignment 10:")
print("\nSequence 1: close-first")
print("Sequence 2: close-second")
get_output(x,y)
"""

"""
# Linear-distant
x = fasta("Assignment1/input files/distant-first.fasta").get("2")
y = fasta("Assignment1/input files/distant-second.fasta").get("2")
print("\nAlignment 1:")
print("\nSequence 1: distant-first")
print("Sequence 2: distant-second")
get_output(x, y)

x = fasta("Assignment1/input files/distant-first.fasta").get("4")
y = fasta("Assignment1/input files/distant-second.fasta").get("4")
print("\nAlignment 2:")
print("\nSequence 1: distant-first")
print("Sequence 2: distant-second")
get_output(x, y)

x = fasta("Assignment1/input files/distant-first.fasta").get("9")
y = fasta("Assignment1/input files/distant-second.fasta").get("9")
print("\nAlignment 3:")
print("\nSequence 1: distant-first")
print("Sequence 2: distant-second")
get_output(x, y)

x = fasta("Assignment1/input files/distant-first.fasta").get("10")
y = fasta("Assignment1/input files/distant-second.fasta").get("10")
print("\nAlignment 4:")
print("\nSequence 1: distant-first")
print("Sequence 2: distant-second")
get_output(x, y)

x = fasta("Assignment1/input files/distant-first.fasta").get("12")
y = fasta("Assignment1/input files/distant-second.fasta").get("12")
print("\nAlignment 5:")
print("\nSequence 1: distant-first")
print("Sequence 2: distant-second")
get_output(x, y)

x = fasta("Assignment1/input files/distant-first.fasta").get("19")
y = fasta("Assignment1/input files/distant-second.fasta").get("19")
print("\nAlignment 6:")
print("\nSequence 1: distant-first")
print("Sequence 2: distant-second")
get_output(x, y)

x = fasta("Assignment1/input files/distant-first.fasta").get("21")
y = fasta("Assignment1/input files/distant-second.fasta").get("21")
print("\nAlignment 7:")
print("\nSequence 1: distant-first")
print("Sequence 2: distant-second")
get_output(x, y)

x = fasta("Assignment1/input files/distant-first.fasta").get("24")
y = fasta("Assignment1/input files/distant-second.fasta").get("24")
print("\nAlignment 8:")
print("\nSequence 1: distant-first")
print("Sequence 2: distant-second")
get_output(x, y)

x = fasta("Assignment1/input files/distant-first.fasta").get("27")
y = fasta("Assignment1/input files/distant-second.fasta").get("27")
print("\nAlignment 9:")
print("\nSequence 1: distant-first")
print("Sequence 2: distant-second")
get_output(x, y)

x = fasta("Assignment1/input files/distant-first.fasta").get("34")
y = fasta("Assignment1/input files/distant-second.fasta").get("34")
print("\nAlignment 10:")
print("\nSequence 1: distant-first")
print("Sequence 2: distant-second")
get_output(x, y)
"""
