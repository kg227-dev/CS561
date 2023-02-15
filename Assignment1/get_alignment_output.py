from linear import *
from affine import *


def sequence_alignment(linex, liney, f, score=0, affine=0):
    if affine == 0:
        linex, liney, score = needleman_wunsch(linex, liney)

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
    print("Matches: " + str(matches), file=f)
    print("Percent identity: " + str(percent_identity)+"%", file=f)
    print("Indels: " + "number="+str(indels) + ", mean length=" +
          str(round(get_mean_indel_length(out_line), 1)), file=f)
    print("Alignment length: " + str(max(len(linex), len(liney))), file=f)
    print("Score="+str(score), file=f)
    for i in range(math.ceil(max(len(linex), len(liney))/60)):
        print("\n"+linex[60*i:60*(i+1)], file=f)
        print(out_line[60*i:60*(i+1)], file=f)
        print(liney[60*i:60*(i+1)], file=f)


def output(file1, file2, save_to_file, matrix, is_affine=0):
    file_output = open(save_to_file, 'x')

    # used for titling in distant text file genes
    distant = ["2", "4", "9", "10", "12", "19", "21", "24", "27", "34"]

    count = 1

    #######
    # Uncommented for different gene titling from close vs distant text file formatting
    #######
    # for _ in range(10):
    # first = "first" + str(count)
    # second = "second" + str(count)
    for gene in distant:
        first = gene
        second = gene

        # Get one sequence from each file
        x = fasta(file1).get(first)
        y = fasta(file2).get(second)
        print(x)
        # Get dictionary of nucletide match/mismatch score
        score_matrix = get_scoring_matrix_dict(matrix)

        # perform DEFAULT linear (is_affine=0) or affine (is_affine=1) alignment
        if is_affine:
            alignment_1, alignment_2, affine_align_score = affine(
                x, y, score_matrix)
        else:
            alignment_1, alignment_2 = needleman_wunsch(x, y)

        print("\nAlignment {}:".format(count), file=file_output)

        #######
        # Uncommented for different gene titling from close vs distant text file formatting
        #######
        # print("\nSequence #1: close-first", file=file_output)
        # print("Sequence #2: close-second", file=file_output)

        print("\nSequence #1: distant-first ({})".format(gene), file=file_output)
        print("Sequence #2: distant-second ({})".format(gene), file=file_output)

        sequence_alignment(alignment_1, alignment_2,
                           file_output, affine_align_score, 1)

        count += 1

    file_output.close

    return
