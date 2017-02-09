# -*- coding: utf-8 -*-

def get_alignments(iseq, jseq, backtracking, end_cell):
    iseq_r = []
    jseq_r = []
    current_cell_v = backtracking[end_cell[0]][end_cell[1]]
    current_cell_index = end_cell
    while current_cell_v != "s":
        if current_cell_v == "u":
            jseq_r.append("-")
            iseq_r.append(iseq[current_cell_index[0] - 1])
            current_cell_index = (current_cell_index[0] - 1, current_cell_index[1])
        elif current_cell_v == "l":
            jseq_r.append("-")
            jseq_r.append(jseq_r[current_cell_index[1] - 1])
            current_cell_index = (current_cell_index[0], current_cell_index[1] - 1)
        else:
            iseq_r.append(iseq[current_cell_index[0] - 1])
            jseq_r.append(jseq[current_cell_index[1] - 1])
            current_cell_index = (current_cell_index[0] - 1, current_cell_index[1] - 1)
        current_cell_v = backtracking[current_cell_index[0]][current_cell_index[1]]
    return "".join(reversed(iseq_r)), "".join(reversed(jseq_r))


def get_consensus_sequences(aligned_seqs):
    result = []
    for fl, sl in zip(*aligned_seqs):
        if fl == "-":
            candidate = sl
        elif sl == "-":
            candidate = fl
        elif fl == "N":
            candidate = sl
        elif sl == "N":
            candidate = fl
        elif sl == fl:
            candidate = sl
        else:
            candidate = "N"
        result.append(candidate)
    return "".join(result)


def substitute_score(letter1, letter2):
    if letter1 == "N" or letter2 == "N" or letter1 == letter2:
        return 0
    return 1


def bounded_alignment(seq1, seq2, max_edit_distance=-1):
    if max_edit_distance == -1:
        max_edit_distance = max(len(seq1), len(seq2))
    i_range = len(seq1) + 1
    j_range = len(seq2) + 1
    M = [[0] * j_range for _ in range(i_range)]
    Mbt = [["s" for _ in range(j_range)] for _ in range(i_range)]
    for j in range(1, j_range):
        M[0][j] = j
        Mbt[0][j] = "l"
    for i in range(1, i_range):
        M[i][0] = i
        Mbt[i][0] = "u"
    for i in range(1, i_range):
        for j in range(1, j_range):
            if abs(i - j) > max_edit_distance:
                continue
            options = [(M[i - 1][j - 1] + substitute_score(seq1[i - 1], seq2[j - 1]), "d")]
            if abs(i - j) < max_edit_distance:
                options.append((M[i - 1][j] + 1, "u"))
                options.append((M[i][j - 1] + 1, "l"))
            else:
                if i > j:
                    options.append((M[i - 1][j] + 1, "u"))
                else:
                    options.append((M[i][j - 1] + 1, "l"))
            result = min(options, key=lambda entry: entry[0])
            M[i][j] = result[0]
            Mbt[i][j] = result[1]
    if M[-1][-1] > max_edit_distance:
        return AlignmentResult(seq1=seq1, seq2=seq2, consensus=None, al_seq1=None, al_seq2=None, edit_distance=M[-1][-1])
    aligned_seqs = get_alignments(iseq=seq1, jseq=seq2, backtracking=Mbt, end_cell=(i_range - 1, j_range - 1))
    consensus_seq = get_consensus_sequences(aligned_seqs=aligned_seqs)
    return AlignmentResult(seq1=seq1, seq2=seq2, consensus=consensus_seq, al_seq1=aligned_seqs[0], al_seq2=aligned_seqs[1], edit_distance=M[-1][-1])


class AlignmentResult(object):
    def __init__(self, seq1, seq2, consensus, al_seq1, al_seq2, edit_distance):
        self.seq1 = seq1
        self.seq2 = seq2
        self.al_seq1 = al_seq1
        self.al_seq2 = al_seq2
        self.consensus = consensus
        self.edit_distance = edit_distance
