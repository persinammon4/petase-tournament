import pandas as pd
from typing import List, Tuple
from Bio import Align

def load_wildtype(filepath: str) -> str:
    """
    Loads the wildtype sequence from the CSV file.
    Assumes the first row contains the WT sequence in 'Wt AA Sequence' column.
    """
    df = pd.read_csv(filepath)
    # Based on file read, the column name is 'Wt AA Sequence'
    if 'Wt AA Sequence' not in df.columns:
        raise ValueError(f"Column 'Wt AA Sequence' not found in {filepath}")
    return df.iloc[0]['Wt AA Sequence']

def load_test_data(filepath: str) -> pd.DataFrame:
    """
    Loads the test dataset.
    """
    return pd.read_csv(filepath)

def get_mutations(wt_seq: str, mutant_seq: str) -> List[str]:
    """
    Identifies mutations between WT and mutant sequences.
    Uses pairwise alignment to handle indels (insertions/deletions).
    Returns a list of mutation strings (e.g., "A123T", "Indel").
    """
    if wt_seq == mutant_seq:
        return ["WildType"]

    # If lengths differ, use alignment to find the indel
    if len(wt_seq) != len(mutant_seq):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        # Use simple scoring for protein alignment
        aligner.match_score = 1.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = -0.5
        
        alignments = aligner.align(wt_seq, mutant_seq)
        # We just need to know it's an indel for logging purposes
        # The scoring logic (PLL) handles the whole sequence, so exact indel parsing isn't critical for the score itself
        return [f"Indel_L{len(wt_seq)}_to_L{len(mutant_seq)}"]

    # Fast path for substitutions (equal length)
    mutations = []
    for i, (wt, mut) in enumerate(zip(wt_seq, mutant_seq)):
        if wt != mut:
            # 1-based indexing for standard notation
            mutations.append(f"{wt}{i+1}{mut}")
            
    return mutations
