import numpy as np
from parsing import *
gap_pen = -1

def get_score(a,b):
    if a == b:
        score = 5
    elif((a=='A') and (b=='G')) or ((a=='G') and (b=='A')):
        score = -1
    elif((a=='C') and (b=='T')) or ((a=='T') and (b=='C')):
        score = -1
    elif((a=='A') and (b=='T')) or ((a=='T') and (b=='A')):
        score = -4
    elif((a=='C') and (b=='G')) or ((a=='G') and (b=='C')):
        score = -4
    elif a=='-' and b=='-':
        score = -1
    else:
        score = -2
    return score

def needleman_wunsch(x, y):
    
    #find length of both sequences
    lenx = len(x)
    leny = len(y)

    #create matrix of zeros
    F = np.zeros((leny+1, lenx+1))

    #first column
    for i in range(0, leny+1):
        F[i][0] = gap_pen * i

    #first row
    for j in range(0, lenx+1):
        F[0][j] = gap_pen * j 

    #fill values in score matrix 
    for i in range(1, len(y)+1):
        for j in range(1, len(x)+1):
            #check top, left, and diag cells
            match = F[i-1,j-1]+ get_score(x[j-1],y[i-1])
            deletion = F[i-1][j] + gap_pen
            insertion = F[i][j-1] + gap_pen
            #get the max score
            F[i,j] = max(match, deletion, insertion)

    print(F)
    #store the alignment 
    linex = ""
    liney = ""

    #set position of i and j 
    i = leny
    j = lenx

    #end traceback when matrix reaches the top or the left side
    while i>0 and j>0: 
        current = F[i][j]
        diag = F[i-1][j-1]
        up = F[i][j-1]
        left = F[i-1][j]
        if current == diag + get_score(x[j-1], y[i-1]):
            linex += x[j-1]
            liney += y[i-1]
            i -=1
            j -=1
        elif current == up + gap_pen:
            linex += x[j-1]
            liney += '-'
            j -=1
        elif current == left + gap_pen:
            linex += '-'
            liney += y[i-1]
            i -= 1
    
    #trace up to the top left cell
    while j>0:
        linex += x[j-1]
        liney += '-'
        j-=1
    while i>0:
        linex += '-'
        liney += y[i-1]
        i-=1
    
    #traverse the score matrix from the bottom right, our two sequence will be reversed
    linex = linex[::-1]
    liney = liney[::-1]

    return(linex, liney)

x = fasta("Assignment1/close-first.fasta").get("first1")[0:40]
y = fasta("Assignment1/close-second.fasta").get("second1")[0:40]


linex, liney = needleman_wunsch(x,y)
print(linex + "\n" + liney)
    