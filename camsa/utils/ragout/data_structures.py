# -*- coding: utf-8 -*-

from camsa.core.data_structures import Sequence


class Block(object):
    def __init__(self, name, start, end, strand, parent_seq=None, annotation=None):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.parent_seq = parent_seq
        self._annotation = annotation

    @property
    def length(self):
        return self.end - self.start

    @property
    def annotation_name(self):
        if self._annotation is None:
            return ".".join([str(self.parent_seq.genome_name), str(self.parent_seq.seq_name), str(self.name), str(self.start) + "-" + str(self.end), str(self.length)])
        return self._annotation

    def __str__(self):
        return self.annotation_name

    def __repr__(self):
        return str(self)


class RagoutSequence(Sequence):
    def __init__(self, name, ragout_id, length=None):
        super(RagoutSequence, self).__init__(name=name, length=length)
        self.ragout_id = ragout_id

    @property
    def genome_name(self):
        return self.name.split(".")[0]

    @property
    def seq_name(self):
        return self.name.split(".", 1)[1]


