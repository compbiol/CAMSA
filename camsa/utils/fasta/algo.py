# -*- coding: utf-8 -*-
import sys

ed_matrix = {
    "A": {"A": 0, "C": 1, "G": 1, "T": 1, "N": 0},
    "C": {"A": 1, "C": 0, "G": 1, "T": 1, "N": 0},
    "G": {"A": 1, "C": 1, "G": 0, "T": 1, "N": 0},
    "T": {"A": 1, "C": 1, "G": 1, "T": 0, "N": 0},
    "N": {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0},
    "gap": 1,
}

al_matrix = {
    "A": {"A": 0, "C": 1, "G": 1, "T": 1, "N": 0},
    "C": {"A": 1, "C": 0, "G": 1, "T": 1, "N": 0},
    "G": {"A": 1, "C": 1, "G": 0, "T": 1, "N": 0},
    "T": {"A": 1, "C": 1, "G": 1, "T": 0, "N": 0},
    "N": {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0},
    "gap": 1,
}


def get_alignments(iseq, jseq, backtracking, end_cell):
    iseq_r = []
    jseq_r = []
    current_cell_v = backtracking[end_cell[0]][end_cell[1]]
    current_cell_index = end_cell
    while current_cell_v != "s":
        if current_cell_v == "u":
            jseq_r.append("-")
            iseq_r.append(iseq[current_cell_index[0]])
            current_cell_index = (current_cell_index[0] - 1, current_cell_index[1])
        elif current_cell_v == "l":
            jseq_r.append("-")
            jseq_r.append(jseq_r[current_cell_index[1]])
            current_cell_index = (current_cell_index[0], current_cell_index[1] - 1)
        else:
            iseq_r.append(iseq[current_cell_index[0]])
            jseq_r.append(jseq[current_cell_index[1]])
            current_cell_index = (current_cell_index[0] - 1, current_cell_index[1] - 1)
    return "".join(reversed(iseq_r)), "".join(reversed(jseq_r))


def get_consensus_sequences(aligned_seqs, first):
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
        else:
            candidate = fl if first else sl
        result.append(candidate)
    return "".join(result)


def sp_alignment(suffix, prefix,
                 ed_matrix=ed_matrix,
                 min_suffix=0, min_prefix=0,
                 max_edit=-1):
    """
    :param prefix:
    :param suffix:
    :param matrix:
    :return: consensus seq, aligned seqs (tuple of two strings read w.r.t. produced alignment), edit distance
    """
    if max_edit == - 1:
        max_edit = sys.maxsize
    i_range = len(suffix) + 1
    j_range = len(prefix) + 1
    M = [[0] * j_range for _ in range(i_range)]
    for i in range(min_suffix):
        M[i][0] = sys.maxsize
    for j in range(min_prefix):
        M[-1][j] = sys.maxsize
    Mbt = [["s" for _ in range(j_range)] for _ in range(i_range)]
    for j in range(j_range):
        M[0][j] = ed_matrix["gap"] * j
        Mbt[0][j] = "l"
    for i in range(1, i_range-1):
        for j in range(j_range-1):
            suffix_letter = suffix[i]
            prefix_letter = prefix[j]
            up = M[i - 1][j] + ed_matrix["gap"]
            diag = M[i][j] + ed_matrix[suffix_letter][prefix_letter]
            left = M[i][j - 1] + ed_matrix["gap"]
            M[i][j] = max(diag, up, left)
            if M[i][j] == diag:
                Mbt[i][j] = "d"
            elif M[i][j] == up:
                Mbt[i][j] = "u"
            else:
                Mbt[i][j] = "l"

    ed_scores = [(i, value) for i, value in enumerate(M[i_range - 1][::-1]) if value <= max_edit]
    if len(ed_scores) == 0:
        return None, None, None
    best_entry = sorted(ed_scores, key=lambda entry: (entry[1]))[-1]
    aligned_seqs = get_alignments(iseq=suffix, jseq=prefix,
                                  backtracking=Mbt,
                                  end_cell=(i_range - 1, best_entry[1]))
    consensus_seq = get_consensus_sequences(aligned_seqs=aligned_seqs, first=True)
    return consensus_seq, aligned_seqs, best_entry[1]


def bounded_alignment(seq1, seq2, matrix):
    pass























