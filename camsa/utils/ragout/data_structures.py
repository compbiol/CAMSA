# -*- coding: utf-8 -*-


class Block(object):
    def __init__(self, name, start, end, strand, parent_seq=None):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.parent_seq = parent_seq


def parse_ragout_block_entry(entry_as_a_string, sequences_by_ids=None):
    pass
