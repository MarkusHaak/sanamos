import pytest
import numpy as np
from sanamos import *

m1, seq1 = "AGCT", "AAAGCTGG"

def test_empty_sa():
    """Minimum length of a sequence is 2!"""
    with pytest.raises(ValueError) as e_info:
        sa = get_suffix_array("empty", "")

@pytest.mark.parametrize('args', 
    [
        # longer motif than seq
        ("AAAAAA", "AAAA", 0, [], []),
        # different motifs
        ("AGCT", "AAAGCTGG", 0, [2], [5]),
        ("GC", "AGCGCT", 0, [1,3], [2,4]),
        ("A", "AAAA", 0, [0,1,2,3], []),
        ("ANA", "AAAA", 0, [0,1], []),
        ("ANA", "AATA", 0, [1], []),
        # change position of interest
        ("AGCT", "AAAGCTGG", 1, [3], [4]),
        ("GC", "AGCGCT", 1, [2,4], [1,3]),
        ("ANA", "AAAA", 1, [1,2], []),
        ("ANA", "AATA", 1, [2], []),
    ])
def test_motif_search(args):
    motif, seq, poi, exp_fwd, exp_rev = args
    sa = get_suffix_array("test", seq, "")
    idx_fwd, idx_rev = find_motif(motif, sa, poi=poi)
    assert (np.array_equal(np.sort(idx_fwd), np.array(exp_fwd)) and 
            np.array_equal(np.sort(idx_rev), np.array(exp_rev)))