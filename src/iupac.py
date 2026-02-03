"""IUPAC ambiguity code handling for PAM sequences."""

from itertools import product
from typing import List

# IUPAC nucleotide ambiguity codes
IUPAC = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'N': ['A', 'C', 'G', 'T'],
    'R': ['A', 'G'],      # puRine
    'Y': ['C', 'T'],      # pYrimidine
    'S': ['G', 'C'],      # Strong
    'W': ['A', 'T'],      # Weak
    'K': ['G', 'T'],      # Keto
    'M': ['A', 'C'],      # aMino
    'B': ['C', 'G', 'T'], # not A
    'D': ['A', 'G', 'T'], # not C
    'H': ['A', 'C', 'T'], # not G
    'V': ['A', 'C', 'G'], # not T
}

# Complement bases
COMPLEMENT = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D',
}


def expand_iupac(pattern: str) -> List[str]:
    """
    Expand an IUPAC pattern to all concrete sequences.

    Args:
        pattern: PAM pattern with IUPAC codes (e.g., "NGG")

    Returns:
        List of all concrete sequences (e.g., ["AGG", "CGG", "GGG", "TGG"])
    """
    pattern = pattern.upper()

    # Get possible bases for each position
    options = []
    for char in pattern:
        if char not in IUPAC:
            raise ValueError(f"Unknown IUPAC code: {char}")
        options.append(IUPAC[char])

    # Generate all combinations
    return [''.join(combo) for combo in product(*options)]


def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence (can include IUPAC codes)

    Returns:
        Reverse complement sequence
    """
    seq = seq.upper()
    complement = []
    for base in seq:
        if base not in COMPLEMENT:
            raise ValueError(f"Unknown base: {base}")
        complement.append(COMPLEMENT[base])
    return ''.join(reversed(complement))


def generate_all_pams(length: int) -> List[str]:
    """
    Generate all possible PAM sequences of a given length.

    Args:
        length: Length of PAM sequences to generate

    Returns:
        List of all PAM sequences (e.g., for length=2: ["AA", "AC", ...])
    """
    bases = ['A', 'C', 'G', 'T']
    return [''.join(combo) for combo in product(bases, repeat=length)]
