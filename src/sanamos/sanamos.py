import os
import random
import sys
import ctypes as ct
from ctypes import POINTER as PT

import numpy as np
from Bio import SeqIO
import pathlib

comp_trans = str.maketrans("ACGTMRWSYKVHDBN", "TGCAKYWSRMBDHVN")

class SuffixArray(ct.Structure):
    """ creates an object to match struct SuffixArray """

    _fields_ = [('sa', PT(ct.c_int)),
                ('s', PT(ct.c_int)),
                ('lcp', PT(ct.c_int)),
                ('rmq', PT(PT(ct.c_int))),
                ('sar', PT(ct.c_int)),
                ('index', PT(PT(PT(ct.c_int)))),
                ('n', ct.c_int),
                ('K', ct.c_int)]

# type for numpy nd-array pointer
nd_pp = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
single_2d_pp = np.ctypeslib.ndpointer(dtype=np.single, ndim=2, flags='C_CONTIGUOUS')
single_1d_pp = np.ctypeslib.ndpointer(dtype=np.single, ndim=1, flags='C_CONTIGUOUS')
char_2d_pp = np.ctypeslib.ndpointer(dtype=np.byte, ndim=2, flags='C_CONTIGUOUS')
int_2d_pp = np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='C_CONTIGUOUS')

# Load the shared library into ctypes
libfile = pathlib.Path(__file__).parent / "sa.so"
c_lib = ct.CDLL(str(libfile))
c_lib.encode_str.restype = ct.c_int
c_lib.encode_str.argtypes = [
    ct.c_char_p, PT(ct.c_int), PT(ct.c_int)]
c_lib.get_na_enc.argtypes = []
c_lib.get_na_enc.restype = PT(ct.c_int)
c_lib.na_enc_init.argtypes = []
c_lib.skew.argtypes = [
    PT(ct.c_int), PT(ct.c_int), ct.c_int, ct.c_int]
c_lib.kasai.argtypes = [
    PT(ct.c_int), PT(ct.c_int), ct.c_int, PT(ct.c_int), PT(ct.c_int)]
c_lib.construct_rmq_array.argtypes = [
    PT(ct.c_int), ct.c_int, PT(PT(ct.c_int))]
c_lib.create_index.argtypes = [
  PT(ct.c_int), PT(ct.c_int), ct.c_int, PT(ct.c_int), PT(PT(PT(ct.c_int)))]
c_lib.find_motif.argtypes = [
    ct.c_char_p, ct.c_int, PT(ct.c_int), 
    PT(ct.c_int), PT(ct.c_int), ct.c_int, 
    PT(PT(ct.c_int)), PT(PT(PT(ct.c_int))), PT(ct.c_int),
    PT(PT(ct.c_int))]
c_lib.motif_means.argtypes = [
    char_2d_pp, ct.c_int, ct.c_int, ct.c_char_p,
    PT(ct.c_void_p), PT(ct.c_void_p),
    PT(PT(SuffixArray)), ct.c_int,
    single_2d_pp, int_2d_pp]
c_lib.motif_medians.argtypes = [
    char_2d_pp, ct.c_int, ct.c_int, ct.c_char_p,
    PT(ct.c_void_p), PT(ct.c_void_p),
    PT(PT(SuffixArray)), ct.c_int,
    single_2d_pp, int_2d_pp]
c_lib.all_motif_medians.argtypes = [
    char_2d_pp, ct.c_int, ct.c_int, ct.c_int,
    single_1d_pp, single_1d_pp,
    PT(ct.c_int), PT(ct.c_int), PT(ct.c_int), ct.c_int, 
    PT(PT(ct.c_int)), PT(PT(PT(ct.c_int))), PT(ct.c_int),
    single_2d_pp, int_2d_pp]
c_lib.get_indices.argtypes = [
    ct.c_int, ct.c_int, PT(ct.c_int), nd_pp]
c_lib.get_indices_from_bipartite_search.argtypes = [
    PT(ct.c_int), ct.c_int, nd_pp]
c_lib.suffix_array_init.argtypes = [
    ct.c_char_p, PT(ct.c_int), PT(SuffixArray)]
c_lib.delete_suffix_array.argtypes = [
    PT(SuffixArray)]
c_lib.print_sa_lcd.argtypes = [
    PT(SuffixArray), ct.c_char_p]
c_lib.write_sa.argtypes = [
    PT(SuffixArray), ct.c_char_p]
c_lib.read_sa.argtypes = [
    ct.c_char_p, ct.c_int, PT(SuffixArray)]

# initialize library
c_lib.na_enc_init()

def get_suffix_array(record_id, seq, cache_fp=''):
    if len(seq) <= 1:
        # actually, the minimum seq length is the index size, which is 8
        raise ValueError("Minimum sequence length for suffix array creation is 1.") 
    sa = SuffixArray()
    if cache_fp:
        if not os.path.exists(os.path.dirname(cache_fp)):
            os.makedirs(os.path.dirname(cache_fp))
        c_cache_fp = ct.c_char_p(cache_fp.encode('utf8'))
        if os.path.exists(cache_fp):
            c_lib.read_sa(c_cache_fp, len(seq) + 1, ct.byref(sa))
            return sa
    target = str(seq).encode('utf8')
    target = ct.c_char_p(target)
    c_lib.suffix_array_init(target, c_lib.get_na_enc(), ct.byref(sa))
    c_lib.skew(sa.s, sa.sa, sa.n, sa.K)
    c_lib.kasai(sa.s, sa.sa, sa.n, sa.lcp, sa.sar)
    c_lib.construct_rmq_array(sa.lcp, sa.n, sa.rmq)
    c_lib.create_index(sa.s, sa.sa, sa.n, sa.lcp, sa.index)
    if cache_fp:
        c_lib.write_sa(ct.byref(sa), c_cache_fp)
    return sa

def delete_suffix_array(sa):
    c_lib.delete_suffix_array(sa)
    del sa

def print_sa_lcd(sa, seq):
    c_str = ct.c_char_p(seq.encode('utf8'))
    c_lib.print_sa_lcd(ct.byref(sa), c_str)

def reverse_complement(seq):
    return seq.translate(comp_trans)[::-1]

def find_motif(motif, sa, poi=0, rc=True):
    if len(motif) >= sa.n:
        return np.array([]), np.array([])
    c_query = ct.c_char_p(motif.upper().encode('utf8'))
    if rc:
        motif_rc = reverse_complement(motif.upper())
        c_query_rc = ct.c_char_p(motif_rc.encode('utf8'))
    c_indices = PT(ct.c_int)()
    n = c_lib.find_motif(c_query, len(motif),
                sa.sa, sa.lcp, sa.s, sa.n, sa.rmq, sa.index, sa.sar,
                ct.byref(c_indices))
    if n:
        # TODO: think about converting c array directly ("by buffer" numpy function)
        # problem: memory is still managed by c then and needs to be freed manually
        indices = np.empty(shape=(n,), dtype=np.int32)
        c_lib.get_indices_from_bipartite_search(c_indices, n, indices)
        if poi:
            indices += poi
    else:
        indices = np.array([])
    if rc:
        if motif_rc == motif.upper():
            # for palindromic motifs, we don't need to search again
            indices_rc = indices.copy()
            indices_rc += len(motif) - 2*poi - 1
        else:
            c_indices_rc = PT(ct.c_int)()
            n = c_lib.find_motif(c_query_rc, len(motif),
                    sa.sa, sa.lcp, sa.s, sa.n, sa.rmq, sa.index, sa.sar,
                    ct.byref(c_indices_rc))
            if n:
                indices_rc = np.empty(shape=(n,), dtype=np.int32)
                c_lib.get_indices_from_bipartite_search(c_indices_rc, n, indices_rc)
                indices_rc += len(motif) - poi - 1
            else:
                indices_rc = np.array([])
        return indices, indices_rc
    else:
        return indices

def motif_means(motifs, max_mlen, fwds, revs, SAs, bases="AC"):
    # construct motif array
    motif_array = np.array([m.ljust(max_mlen+1, '\0').encode('utf8') for m in motifs], dtype=f'|S{max_mlen+1}')
    motif_array = motif_array.view(np.byte).reshape((motif_array.size, -1))
    means = np.full((len(motifs), max_mlen), dtype=np.single, fill_value=np.nan)
    counts = np.full((len(motifs), max_mlen), dtype=np.int32, fill_value=0)
    c_bases = ct.c_char_p((bases + '\0').encode('utf8'))

    fwd_ = [np.array(fwd, dtype=np.single) for fwd in fwds]
    rev_ = [np.array(rev, dtype=np.single) for rev in revs]
    fwd_arr = (ct.c_void_p * len(SAs))(*[ct.cast(np.ctypeslib.as_ctypes(fwd), ct.c_void_p) for fwd in fwd_])
    rev_arr = (ct.c_void_p * len(SAs))(*[ct.cast(np.ctypeslib.as_ctypes(rev), ct.c_void_p) for rev in rev_])
    SAs_arr = (PT(SuffixArray) * len(SAs))(*[ct.cast(ct.byref(sa), PT(SuffixArray)) for sa in SAs])
    c_lib.motif_means(
        motif_array, len(motifs), max_mlen+1, c_bases,
        fwd_arr, rev_arr,
        SAs_arr, len(SAs),
        means, counts)
    return means, counts

def motif_medians(motifs, max_mlen, fwds, revs, SAs, bases="AC"):
    if type(SAs) != list:
        fwds = [fwds]
        revs = [revs]
        SAs = [SAs]
    # construct motif array
    motif_array = np.array([m.ljust(max_mlen+1, '\0').encode('utf8') for m in motifs], dtype=f'|S{max_mlen+1}')
    motif_array = motif_array.view(np.byte).reshape((motif_array.size, -1))
    medians = np.full((len(motifs), max_mlen), dtype=np.single, fill_value=np.nan)
    counts = np.full((len(motifs), max_mlen), dtype=np.int32, fill_value=0)
    c_bases = ct.c_char_p((bases + '\0').encode('utf8'))
    
    fwd_ = [np.array(fwd, dtype=np.single) for fwd in fwds]
    rev_ = [np.array(rev, dtype=np.single) for rev in revs]
    fwd_arr = (ct.c_void_p * len(SAs))(*[ct.cast(np.ctypeslib.as_ctypes(fwd), ct.c_void_p) for fwd in fwd_])
    rev_arr = (ct.c_void_p * len(SAs))(*[ct.cast(np.ctypeslib.as_ctypes(rev), ct.c_void_p) for rev in rev_])
    SAs_arr = (PT(SuffixArray) * len(SAs))(*[ct.cast(ct.byref(sa), PT(SuffixArray)) for sa in SAs])
    c_lib.motif_medians(
        motif_array, len(motifs), max_mlen+1, c_bases,
        fwd_arr, rev_arr,
        SAs_arr, len(SAs),
        medians, counts)
    return medians, counts

def all_motif_medians(motifs, max_mlen, fwd, rev, sa, pad=2):
    # construct motif array
    motif_array = np.array([m.ljust(max_mlen+1, '\0').encode('utf8') for m in motifs], dtype=f'|S{max_mlen+1}')
    motif_array = motif_array.view(np.byte).reshape((motif_array.size, -1))
    medians = np.full((len(motifs), max_mlen+2*pad), dtype=np.single, fill_value=np.nan)
    counts = np.full((len(motifs), max_mlen+2*pad), dtype=np.int32, fill_value=0)
    fwd_ = np.array(fwd, dtype=np.single)
    rev_ = np.array(rev, dtype=np.single)
    c_lib.all_motif_medians(
        motif_array, len(motifs), max_mlen+1, pad,
        fwd_, rev_,
        sa.sa, sa.lcp, sa.s, sa.n, 
        sa.rmq, sa.index, sa.sar,
        medians, counts)
    return medians, counts
