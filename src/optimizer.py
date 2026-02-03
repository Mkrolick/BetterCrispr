"""Find optimal PAM sequences to maximize genome coverage."""

from typing import List, Tuple, Dict, Optional
from tqdm import tqdm

from .iupac import generate_all_pams, expand_iupac
from .pam_finder import (find_pam_sites, get_target_regions,
                         COVERAGE_CUT_SITE, COVERAGE_EDITING_WINDOW, COVERAGE_GUIDE_BINDING)
from .coverage import CoverageTracker


def find_best_next_pam(
    genome_sequences: Dict[str, str],
    current_coverage: CoverageTracker,
    pam_length: int = 3,
    top_n: int = 10,
    coverage_mode: str = COVERAGE_EDITING_WINDOW,
    show_progress: bool = True
) -> List[Tuple[str, int, float]]:
    """
    Find the best next PAM sequence to maximize additional coverage.

    Args:
        genome_sequences: Dictionary of {chromosome: sequence}
        current_coverage: Current coverage tracker
        pam_length: Length of PAM to search for
        top_n: Number of top PAMs to return
        coverage_mode: Coverage calculation mode ('cut_site', 'editing', or 'guide')
        show_progress: Whether to show progress bar

    Returns:
        List of (pam, additional_bp, additional_pct) tuples, sorted by additional coverage
    """
    candidates = generate_all_pams(pam_length)
    results = []

    total_length = current_coverage.get_total_length()

    iterator = tqdm(candidates, desc="Evaluating PAMs", disable=not show_progress)

    for pam in iterator:
        additional_bp = 0

        for chrom, sequence in genome_sequences.items():
            # Find PAM sites
            sites = find_pam_sites(sequence, pam, search_both_strands=True)

            # Get target regions based on coverage mode
            regions = get_target_regions(sites, len(pam),
                                         coverage_mode=coverage_mode,
                                         seq_length=len(sequence))

            # Count new coverage
            additional_bp += current_coverage.count_new_coverage(chrom, regions)

        additional_pct = (additional_bp / total_length * 100) if total_length > 0 else 0
        results.append((pam, additional_bp, additional_pct))

    # Sort by additional coverage (descending)
    results.sort(key=lambda x: -x[1])

    return results[:top_n]


def analyze_pam_coverage(
    genome_sequences: Dict[str, str],
    pam: str,
    coverage: Optional[CoverageTracker] = None,
    coverage_mode: str = COVERAGE_EDITING_WINDOW,
    show_progress: bool = True
) -> Tuple[int, CoverageTracker]:
    """
    Analyze coverage for a single PAM and update coverage tracker.

    Args:
        genome_sequences: Dictionary of {chromosome: sequence}
        pam: PAM pattern to analyze
        coverage: Existing coverage tracker to update (creates new if None)
        coverage_mode: Coverage calculation mode ('cut_site', 'editing', or 'guide')
        show_progress: Whether to show progress

    Returns:
        Tuple of (total_pam_sites, updated_coverage_tracker)
    """
    if coverage is None:
        coverage = CoverageTracker()
        for chrom, seq in genome_sequences.items():
            coverage.add_chromosome(chrom, len(seq))

    pam_length = len(expand_iupac(pam)[0])  # Get concrete PAM length
    total_sites = 0

    chromosomes = list(genome_sequences.items())
    iterator = tqdm(chromosomes, desc=f"Analyzing {pam}", disable=not show_progress)

    for chrom, sequence in iterator:
        # Find PAM sites
        sites = find_pam_sites(sequence, pam, search_both_strands=True)
        total_sites += len(sites)

        # Get target regions based on coverage mode
        regions = get_target_regions(sites, pam_length,
                                     coverage_mode=coverage_mode,
                                     seq_length=len(sequence))

        # Mark coverage
        coverage.mark_regions_batch(chrom, regions)

    return total_sites, coverage


def analyze_multiple_pams(
    genome_sequences: Dict[str, str],
    pams: List[str],
    coverage_mode: str = COVERAGE_EDITING_WINDOW,
    show_progress: bool = True
) -> Tuple[Dict[str, int], CoverageTracker]:
    """
    Analyze coverage for multiple PAMs.

    Args:
        genome_sequences: Dictionary of {chromosome: sequence}
        pams: List of PAM patterns to analyze
        coverage_mode: Coverage calculation mode ('cut_site', 'editing', or 'guide')
        show_progress: Whether to show progress

    Returns:
        Tuple of (dict mapping PAM to site count, combined coverage tracker)
    """
    coverage = CoverageTracker()
    for chrom, seq in genome_sequences.items():
        coverage.add_chromosome(chrom, len(seq))

    pam_sites = {}

    for pam in pams:
        sites, coverage = analyze_pam_coverage(
            genome_sequences, pam, coverage, coverage_mode, show_progress
        )
        pam_sites[pam] = sites

    return pam_sites, coverage
