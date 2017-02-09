# -*- coding: utf-8 -*-
from Bio.Seq import Seq

from camsa.utils.fasta.algo import bounded_alignment


class GapFilling(object):
    def __init__(self, seq):
        self.seq = seq
        self.score = None

    def suitable_to_fill_the_gap(self, gap_size, gap_size_error):
        return False

    def compute_score(self, gap_size, gap_size_error):
        non_N_content = len([letter for letter in self.seq if letter != "N"])
        self.score = non_N_content


class IntraGapFilling(GapFilling):
    def suitable_to_fill_the_gap(self, gap_size, gap_size_error):
        if gap_size == "?":
            return True
        min_gap_size = gap_size - gap_size_error
        max_gap_size = gap_size + gap_size_error
        return min_gap_size <= len(self.seq) <= max_gap_size


class FlankingGapFilling(GapFilling):
    def __init__(self, seq1, seq2, seq):
        super(FlankingGapFilling, self).__init__(seq)
        self.start_seq = seq1
        self.end_seq = seq2
        self.alignment_result = None

    def suitable_to_fill_the_gap(self, gap_size, gap_size_error, fill_remainder=True):
        if gap_size == "?":
            return False
        min_gap_size = gap_size - gap_size_error
        max_gap_size = gap_size + gap_size_error
        cumulative_length = len(self.start_seq) + len(self.end_seq)
        if cumulative_length < min_gap_size:
            return False
        if cumulative_length > gap_size + 2 * gap_size_error:
            return False
        if cumulative_length > max_gap_size:
            min_overlap = cumulative_length - gap_size - gap_size_error
            max_overlap = cumulative_length - gap_size + gap_size_error
            max_overlap = min(max_overlap, len(self.start_seq), len(self.end_seq))
            l1s = self.start_seq[-max_overlap:]
            l2p = self.end_seq[:max_overlap]
            alignment_result = bounded_alignment(seq1=l1s, seq2=l2p, max_edit_distance=max_overlap - min_overlap)
            if alignment_result.consensus is None:
                return False
            else:
                self.alignment_result = alignment_result
                return True
        if cumulative_length < max_gap_size:
            # there can be an overlap
            return True

    def prepare_seq(self, gap_size, gap_size_error):
        # only called if the gap can be suitable for filling w.r.t. gap size
        self.seq = self.start_seq
        cumulative_length = len(self.start_seq) + len(self.end_seq)
        if cumulative_length < gap_size + gap_size_error:
            self.seq = self.start_seq + self.end_seq
        else:
            max_overlap = cumulative_length - gap_size + gap_size_error
            max_overlap = min(max_overlap, len(self.start_seq), len(self.end_seq))
            self.seq = self.start_seq[:max_overlap] + self.alignment_result.consensus + self.end_seq[max_overlap:]
