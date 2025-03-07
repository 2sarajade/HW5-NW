# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        m, n = len(seqA), len(seqB)
        self._align_matrix = np.full((m+1, n+1), -np.inf)
        self._gapA_matrix = np.full((m+1, n+1), -np.inf)
        self._gapB_matrix = np.full((m+1, n+1), -np.inf)
        self._back = np.full((m+1, n+1), -np.inf)
        self._back_A = np.full((m+1, n+1), -np.inf)
        self._back_B = np.full((m+1, n+1), -np.inf)

        # initialize values
        self._align_matrix[0,0] = 0
        self._back[0,0] = -1
        self._back_A[0,0] = -1
        self._back_B[0,0] = -1
        for i in range(m + 1):
            self._gapA_matrix[i,0] = self.gap_open + (i * self.gap_extend)
            self._back_A[i,0] = 1 
        for j in range(n + 1):
            self._gapB_matrix[0,j] = self.gap_open + (j * self.gap_extend)
            self._back_B[0,j] = 2

        
        # TODO: Implement global alignment here
        for i in range(1, n + 1):
            for j in range(1,m + 1):
                # pick M
                match_val = self.sub_dict.get((seqA[j-1],seqB[i-1]))
                M_match = self._align_matrix[j-1, i-1] + match_val
                M_agap = self._gapA_matrix[j-1, i-1] + match_val
                M_bgap = self._gapB_matrix[j-1, i-1] + match_val #problem here, match val not getting added to bgap value
                M_options = [M_match, M_agap, M_bgap]
                M_choice = np.argmax(M_options)
                # 0 = match, 1 = agap, 2 = bgap, (in event of a tie, match > agap > bgap)
                self._back[j,i ] = M_choice
                self._align_matrix[j,i] = M_options[M_choice] 
                # pick a gap
                A_open = self._align_matrix[j, i-1] + self.gap_open + self.gap_extend
                A_extend = self._gapA_matrix[j, i-1] + self.gap_extend
                A_options = [A_open, A_extend, -np.inf] #same order
                A_choice = np.argmax(A_options)
                self._back_A[j,i] = A_choice
                self._gapA_matrix[j,i] = A_options[A_choice]
                # pick b gap
                B_open = self._align_matrix[j-1, i] + self.gap_open + self.gap_extend
                B_extend = self._gapB_matrix[j-1, i] + self.gap_extend
                B_options = [B_open, -np.inf, B_extend] #same order
                B_choice = np.argmax(B_options)
                self._back_B[j, i] = B_choice
                self._gapB_matrix[j, i] = B_options[B_choice]

        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
 
        m, n = len(self._seqA), len(self._seqB)
        
        prev_direction = 0
        current_mat = self._align_matrix
        current_back = self._back
        # find the start point for backtrace
        if self._gapA_matrix[m,n] >= self._align_matrix[m,n] and self._gapA_matrix[m,n] >= self._gapB_matrix[m,n]:
            prev_direction = 1
            current_back = self._back_A
            current_mat = self._gapA_matrix
            self._back_A
        elif self._gapB_matrix[m,n] >= self._align_matrix[m,n] and self._gapB_matrix[m,n] >= self._gapA_matrix[m,n]:
            prev_direction = 2
            current_back = self._back_B
            current_mat = self._gapB_matrix
            #what I was thinking: I jump the gun on switching matrices (why the q was too early)
            # so setting up in the beginning then the loop will do all the backtrace, not one step done early.
        self.alignment_score = current_mat[m,n]

        while m > 0 and n > 0:
            if prev_direction == 0: #match
                self.seqA_align = self.seqA_align + self._seqA[m-1]
                self.seqB_align = self.seqB_align + self._seqB[n-1]
                current_back = self._back
                prev_direction = current_back[m,n]
                m -= 1
                n -= 1
            elif prev_direction == 1: #a gap
                self.seqA_align = self.seqA_align + "-"
                self.seqB_align = self.seqB_align + self._seqB[n-1]
                current_back = self._back_A
                prev_direction = current_back[m,n]
                n -= 1
            elif prev_direction == 2: # b gap
                self.seqA_align = self.seqA_align + self._seqA[m-1]
                self.seqB_align = self.seqB_align + "-"
                current_back = self._back_B
                prev_direction = current_back[m,n]
                m -= 1
            if prev_direction == -np.inf: #should never happen, here as a debugging tool
                print(m,n)
                print(self._back)
                print(self._back_A)
                print(self._back_B)
                print("********")
                print(self._align_matrix)
                print(self._gapA_matrix)
                print(self._gapB_matrix)
                print("********")
                print(self.seqA_align)
                print(self.seqB_align)
                raise ValueError('-inf pointer')
            
        #deal with gaps at beginning or end
        if n > 0:
            self.seqA_align = self.seqA_align + "-" * n
            self.seqB_align = self.seqB_align + self._seqB[0:n][::-1]
        if m > 0:
            self.seqB_align = self.seqB_align + "-" * m
            self.seqA_align = self.seqA_align + self._seqA[0:m][::-1]
            
            
        #reverse the alignments
        self.seqA_align = self.seqA_align[::-1]
        self.seqB_align = self.seqB_align[::-1]
        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
