"""Find PAM occurrences in genome sequences."""

from typing import List, Tuple, Set
from .iupac import expand_iupac, reverse_complement


# Coverage mode constants
COVERAGE_CUT_SITE = 'cut_site'           # Just the cut position (1bp)
COVERAGE_EDITING_WINDOW = 'editing'       # ~10bp around cut site (practical HDR)
COVERAGE_GUIDE_BINDING = 'guide'          # 20bp guide binding region (original model)
COVERAGE_PRIME_EDITING = 'prime'          # Prime editing window (downstream of nick)

COVERAGE_MODES = {
    COVERAGE_CUT_SITE: {
        'description': 'Only the cut site position (1bp per PAM)',
        'upstream_of_cut': 0,
        'downstream_of_cut': 0,
    },
    COVERAGE_EDITING_WINDOW: {
        'description': 'Efficient HDR editing window (~21bp around cut site)',
        'upstream_of_cut': 10,
        'downstream_of_cut': 10,
    },
    COVERAGE_GUIDE_BINDING: {
        'description': 'Full 20bp guide binding region (not actual editing coverage)',
        'upstream_of_cut': 17,  # Cut is at -3 from PAM, guide extends to -20, so 17bp upstream of cut
        'downstream_of_cut': 2,  # Cut at -3, PAM starts at 0, so 2bp downstream to reach -1
    },
    COVERAGE_PRIME_EDITING: {
        'description': 'Prime editing window (~30bp downstream of nick)',
        'upstream_of_cut': 0,   # Prime editing works downstream of the nick
        'downstream_of_cut': 30, # Typical RT template can cover ~30bp efficiently
    },
}


def find_pam_sites(sequence: str, pam: str, search_both_strands: bool = True) -> List[Tuple[int, str]]:
    """
    Find all PAM sites in a sequence.

    Args:
        sequence: DNA sequence to search
        pam: PAM pattern (can include IUPAC codes)
        search_both_strands: If True, also search reverse complement

    Returns:
        List of (position, strand) tuples where position is the start of the PAM
        strand is '+' for forward, '-' for reverse
    """
    sequence = sequence.upper()
    sites = []

    # Expand PAM to concrete sequences
    concrete_pams = expand_iupac(pam)

    # Search forward strand
    for concrete_pam in concrete_pams:
        pam_len = len(concrete_pam)
        start = 0
        while True:
            pos = sequence.find(concrete_pam, start)
            if pos == -1:
                break
            sites.append((pos, '+'))
            start = pos + 1

    # Search reverse strand
    if search_both_strands:
        # For reverse strand, we search for the reverse complement of the PAM
        rc_pam = reverse_complement(pam)
        concrete_rc_pams = expand_iupac(rc_pam)

        for concrete_pam in concrete_rc_pams:
            pam_len = len(concrete_pam)
            start = 0
            while True:
                pos = sequence.find(concrete_pam, start)
                if pos == -1:
                    break
                sites.append((pos, '-'))
                start = pos + 1

    return sites


def get_cut_sites(pam_sites: List[Tuple[int, str]], pam_length: int,
                  seq_length: int = None) -> List[int]:
    """
    Get the cut site positions for each PAM site.

    Cas9 cuts 3bp upstream of the PAM on the target strand.

    For forward strand (+): cut at pos - 3
    For reverse strand (-): cut at pos + pam_length + 2

    Args:
        pam_sites: List of (position, strand) tuples
        pam_length: Length of the PAM sequence
        seq_length: Total sequence length (for boundary checking)

    Returns:
        List of cut site positions
    """
    cut_sites = []

    for pos, strand in pam_sites:
        if strand == '+':
            # Cut is 3bp upstream of PAM start
            cut_pos = pos - 3
        else:
            # For reverse strand PAM, cut is 3bp after PAM end
            cut_pos = pos + pam_length + 2

        # Boundary checks
        if cut_pos >= 0 and (seq_length is None or cut_pos < seq_length):
            cut_sites.append(cut_pos)

    return cut_sites


def get_target_regions(pam_sites: List[Tuple[int, str]], pam_length: int,
                       coverage_mode: str = COVERAGE_EDITING_WINDOW,
                       seq_length: int = None) -> List[Tuple[int, int]]:
    """
    Get the targetable regions for each PAM site based on coverage mode.

    Coverage modes:
    - 'cut_site': Just the cut position (1bp) - most accurate for "can I edit here?"
    - 'editing': ~10bp around cut site - practical HDR editing window (default)
    - 'guide': 20bp guide binding region - original model (overstates actual coverage)

    For Cas9, the cut site is 3bp upstream of the PAM.

    Args:
        pam_sites: List of (position, strand) tuples
        pam_length: Length of the PAM sequence
        coverage_mode: One of 'cut_site', 'editing', or 'guide'
        seq_length: Total sequence length (for boundary checking)

    Returns:
        List of (start, end) tuples representing target regions (inclusive)
    """
    if coverage_mode not in COVERAGE_MODES:
        raise ValueError(f"Unknown coverage mode: {coverage_mode}. "
                         f"Choose from: {list(COVERAGE_MODES.keys())}")

    mode_config = COVERAGE_MODES[coverage_mode]
    upstream = mode_config['upstream_of_cut']
    downstream = mode_config['downstream_of_cut']

    regions = []

    for pos, strand in pam_sites:
        # First, find the cut site
        if strand == '+':
            cut_pos = pos - 3
        else:
            cut_pos = pos + pam_length + 2

        # Then define region around cut site
        start = cut_pos - upstream
        end = cut_pos + downstream

        # Boundary checks
        if start < 0:
            start = 0
        if seq_length is not None and end >= seq_length:
            end = seq_length - 1

        # Only include valid regions
        if start <= end:
            regions.append((start, end))

    return regions


def get_target_regions_legacy(pam_sites: List[Tuple[int, str]], pam_length: int,
                              target_length: int = 20, seq_length: int = None) -> List[Tuple[int, int]]:
    """
    Legacy function: Get the target regions (20bp upstream of PAM) for each PAM site.

    DEPRECATED: Use get_target_regions with coverage_mode='guide' instead.
    This function is kept for backwards compatibility.

    For forward strand (+): target is [pos - target_length, pos - 1]
    For reverse strand (-): target is [pos + pam_length, pos + pam_length + target_length - 1]

    Args:
        pam_sites: List of (position, strand) tuples
        pam_length: Length of the PAM sequence
        target_length: Length of target region (default 20bp)
        seq_length: Total sequence length (for boundary checking)

    Returns:
        List of (start, end) tuples representing target regions (inclusive)
    """
    regions = []

    for pos, strand in pam_sites:
        if strand == '+':
            # Target is upstream (5') of PAM on forward strand
            start = pos - target_length
            end = pos - 1
        else:
            # Target is downstream of PAM on reverse strand
            # (which is upstream on the complementary strand)
            start = pos + pam_length
            end = pos + pam_length + target_length - 1

        # Boundary checks
        if start < 0:
            start = 0
        if seq_length is not None and end >= seq_length:
            end = seq_length - 1

        # Only include valid regions
        if start <= end:
            regions.append((start, end))

    return regions


def count_pam_sites(sequence: str, pam: str, search_both_strands: bool = True) -> int:
    """
    Count total PAM sites in a sequence.

    Args:
        sequence: DNA sequence to search
        pam: PAM pattern (can include IUPAC codes)
        search_both_strands: If True, also search reverse complement

    Returns:
        Total count of PAM sites
    """
    return len(find_pam_sites(sequence, pam, search_both_strands))
