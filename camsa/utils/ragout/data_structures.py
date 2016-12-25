# -*- coding: utf-8 -*-

from camsa.core.data_structures import Sequence


class Block(object):
    def __init__(self, name, start, end, strand, parent_seq=None):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.parent_seq = parent_seq


class RagoutSequence(Sequence):
    def __init__(self, name, ragout_id, length=None):
        super(RagoutSequence, self).__init__(name=name, length=length)
        self.ragout_id = ragout_id

    @property
    def genome_name(self):
        return self.name.split(".")[0]


