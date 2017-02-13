# -*- coding: utf-8 -*-

from hypothesis import given, settings
from hypothesis.strategies import integers, text, data

from camsa.core.data_structures import Sequence


@settings(perform_health_check=False)
@given(text=text(alphabet="ACGT", min_size=1000, max_size=10000),
       blocks_cnt=integers(min_value=1, max_value=1000),
       data=data())
def test_sequence_reconstruction(text, blocks_cnt, data):
    blocks_lengths = []
    cum_blocks_length = 0
    for i in range(blocks_cnt):
        remaining_length = len(text) - cum_blocks_length - blocks_cnt + i
        length_strategy = integers(min_value=1, max_value=remaining_length)
        block_length = data.draw(length_strategy)
        cum_blocks_length += block_length
        blocks_lengths.append(block_length)
    gap_lengths = []
    cum_gap_length = 0
    for i in range(blocks_cnt):
        remaining_length = len(text) - cum_blocks_length - cum_gap_length
        length_strategy = integers(min_value=0, max_value=remaining_length)
        gap_length = data.draw(length_strategy)
        cum_gap_length += gap_length
        gap_lengths.append(gap_length)
    assert(sum(blocks_lengths) + sum(gap_lengths) <= len(text))
    breakage_cnt = integers(min_value=0, max_value=blocks_cnt-1)
    current_position = 0
    sequences = []
    for cnt in range(blocks_cnt * 2):
        if cnt % 2 == 0:
            # gap
            current_position += gap_lengths[cnt//2]
        else:
            #block
            sequence = Sequence(name=cnt, length=blocks_lengths[cnt//2], parent_seq_id=None,
                                start=current_position, end=current_position+blocks_lengths[cnt//2], strand="+")



