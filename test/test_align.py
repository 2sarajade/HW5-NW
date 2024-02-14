# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    score, a, b = nw.align(seq1, seq2)
    print(nw._align_matrix)
    print(nw._gapA_matrix)
    print(nw._gapB_matrix)
    print("***")
    print(nw._back)
    print(nw._back_A)
    print(nw._back_B)
    assert a == "MYQR"
    assert b == "M-QR" #like this because YQ mismatch value is the same as the gap penalty, and tiebreak has gapa>match>gapb
    
    #assert np.array_equal(nw._align_matrix == [[0, -np.inf, -np.inf, -np.inf, -np.inf],
    #                                           [-np.inf, 5, -12, -12, -14],
    #                                           [-np.inf, -11, 4, 3, -2],
    #                                           [-np.inf, -13, -8, 5, 8]] )
    #print(nw._gapA_matrix)
    #assert np.array_equal(nw._gapA_matrix, [[0, -np.inf, -np.inf, -np.inf, -np.inf],
    #                                           [-11, -12, -2, -3, -4],
    #                                           [-12, -13, -14, -7, -8],
    #                                           [-13, -14, -15, -16, -6]] )
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    score, a, b = nw.align(seq3, seq4)
    #print(nw._align_matrix)
    #assert score == 17
    #assert a == "MAVHQLIRRP"
    #assert b == "M---QLIRHP"




