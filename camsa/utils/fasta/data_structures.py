# -*- coding: utf-8 -*-
from Bio.Seq import Seq


class GapFilling(object):
    def __init__(self, seq):
        self.seq = seq
        self.score = None

    def suitable_to_fill_the_gap(self, gap_size, error_per, error_bp):
        return False

    def compute_score(self, gap_size, error_per, error_bp):
        if gap_size == "?":
            self.score = len(self.seq)
        else:
            self.score = abs(gap_size - len(self.seq))


class IntraGapFilling(GapFilling):
    def suitable_to_fill_the_gap(self, gap_size, error_per, error_bp):
        if gap_size == "?":
            return True
        return gap_size * (1 - error_per / 100) <= len(self.seq) <= gap_size * (1 + error_per / 100) or gap_size - error_bp <= len(self.seq) <= gap_size + error_bp


class FlankingGapFilling(GapFilling):
    def __init__(self, seq1, seq2, seq):
        super(FlankingGapFilling, self).__init__(seq)
        self.start_seq = seq1
        self.end_seq = seq2

    def suitable_to_fill_the_gap(self, gap_size, error_per, error_bp, fill_remainder=True):
        if gap_size == "?":
            return True
        cumulative_length = len(self.start_seq) + len(self.end_seq)
        if cumulative_length < gap_size * (1 - error_per / 100) or cumulative_length < gap_size - error_bp:
            return fill_remainder
        return gap_size * (1 - error_per / 100) <= len(self.seq) <= gap_size * (1 + error_per / 100) or gap_size - error_bp <= len(self.seq) <= gap_size + error_bp

    def prepare_seq(self, gap_size, error_per, error_bp, fill_remainder=True, sep="N"):
        self.seq = self.start_seq
        if gap_size != "?":
            cumulative_length = len(self.start_seq) + len(self.end_seq)
            if cumulative_length < gap_size * (1 - error_per / 100) or cumulative_length < gap_size - error_bp and fill_remainder:
                difference = min(gap_size * (1 - error_per / 100) - cumulative_length, gap_size - error_bp - cumulative_length)
                self.seq += Seq(sep * int(difference))
        self.seq += self.end_seq
