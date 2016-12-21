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


def parse_ragout_seq_entry(entry_as_a_string, seq_id_column=0, seq_length_column=1, seq_description_column=2, delimiter="\t"):
    """
        Seq_id	Size	Description
        27	11234	Anc0.Anc0refChr7043
    """
    data = entry_as_a_string.split(delimiter)
    seq_id = int(float(data[seq_id_column]))
    seq_length = int(float(data[seq_length_column]))
    seq_name = data[seq_description_column]
    return RagoutSequence(name=seq_name, ragout_id=seq_id, length=seq_length)


def parse_ragout_block_entry(entry_as_a_string, sequences_by_ids=None):
    pass


def parse_block_id_string(data_string):
    pass
