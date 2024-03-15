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
        # motifs that are not contained
        ("C", "AAAA", 0, [], []),
        ("GG", "AAAA", 0, [], []),
        # contiguous motifs
        ("AGCT", "AAAGCTGG", 0, [2], [5]),
        ("GC", "AGCGCT", 0, [1,3], [2,4]),
        ("A", "AAAA", 0, [0,1,2,3], []),
        # change position of interest
        ("AGCT", "AAAGCTGG", 1, [3], [4]),
        ("GC", "AGCGCT", 1, [2,4], [1,3]),
        # bipartite motifs
        ("CSNNNNNNNNSCC", "TCCGTTTTTTTGCCCAA", 0, [1,2], []),
        ("CCWGG", "TACCTGGATA", 0, [2], [6]),
        # too long for indexing
        ("AGCTAGCTAGCT", "AGCTAGCTAGCTAGCTAGCT", 0, [0,4,8], [11,15,19]),
        ("AGCTAGCTAGCTNNNNAGCTAGCTTTTT", "AAAAAGCTAGCTAGCTAGCTAGCTAGCTTTTT", 0, [4], [27]),
    ])
def test_motif_search(args):
    motif, seq, poi, exp_fwd, exp_rev = args
    sa = get_suffix_array("test", seq, "")
    idx_fwd, idx_rev = find_motif(motif, sa, poi=poi)
    assert (np.array_equal(np.sort(idx_fwd), np.array(exp_fwd)) and 
            np.array_equal(np.sort(idx_rev), np.array(exp_rev)))

@pytest.mark.parametrize('args', 
    [
        # standard case
        ("GATC", "AAGATCTTGATCAA", 1, np.array([0.,0.,1.,5.,2.,0.,0.,0.,.2,3.,1.,0.,0.,0.]), 
                                      np.array([0.,0.,0.,2.,5.,1.,0.,0.,0.,.2,1.,.1,0.,0.]), 
         np.mean([5., 3., 5., 1.]), 4),
        # nan values
        ("GATC", "AAGATCTTGATCAA", 1, np.array([0.,0.,1.,5.,    2.,0.,0.,0.,.2,3.,1.,0.,0.,0.]), 
                                      np.array([0.,0.,0.,2.,np.nan,1.,0.,0.,0.,.2,1.,.1,0.,0.]), 
         np.mean([5., 3., 1.]), 3),
    ])
def test_motif_means(args):
    motif, seq, poi, fwd, rev, mean, count = args
    sa = get_suffix_array("test", seq, "")
    means, counts = motif_means([motif], len(motif), [fwd], [rev], [sa], bases="AC")
    assert (counts[0][poi] == count and means[0][poi] == mean)

@pytest.mark.parametrize('args', 
    [
        # standard case
        ("GATC", "AAGATCTTGATCAA", 1, np.array([0.,0.,1.,5.,2.,0.,0.,0.,.2,3.,1.,0.,0.,0.]), 
                                      np.array([0.,0.,0.,2.,5.,1.,0.,0.,0.,.2,1.,.1,0.,0.]), 
         np.median([5., 3., 5., 1.]), 4),
         # nan values
        ("GATC", "AAGATCTTGATCAA", 1, np.array([0.,0.,1.,5.,    2.,0.,0.,0.,.2,3.,1.,0.,0.,0.]), 
                                      np.array([0.,0.,0.,2.,np.nan,1.,0.,0.,0.,.2,1.,.1,0.,0.]), 
         np.median([5., 3., 1.]), 3),
    ])
def test_motif_medians(args):
    motif, seq, poi, fwd, rev, mean, count = args
    sa = get_suffix_array("test", seq, "")
    means, counts = motif_medians([motif], len(motif), [fwd], [rev], [sa], bases="AC")
    assert (counts[0][poi] == count and means[0][poi] == mean)
