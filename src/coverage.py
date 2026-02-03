"""Bit-vector coverage tracking for genome positions."""

import numpy as np
from typing import Dict, List, Tuple, Optional


class CoverageTracker:
    """
    Track genome coverage using bit vectors.

    Uses numpy arrays with uint8 dtype for efficient bit operations.
    Each byte stores 8 positions.
    """

    def __init__(self):
        """Initialize empty coverage tracker."""
        self.chromosomes: Dict[str, np.ndarray] = {}
        self.chrom_lengths: Dict[str, int] = {}

    def add_chromosome(self, name: str, length: int):
        """
        Add a chromosome to track.

        Args:
            name: Chromosome name
            length: Length of chromosome in base pairs
        """
        # Allocate bit vector (1 byte per 8 positions)
        num_bytes = (length + 7) // 8
        self.chromosomes[name] = np.zeros(num_bytes, dtype=np.uint8)
        self.chrom_lengths[name] = length

    def mark_position(self, chrom: str, pos: int):
        """Mark a single position as covered."""
        if chrom not in self.chromosomes:
            return

        if pos < 0 or pos >= self.chrom_lengths[chrom]:
            return

        byte_idx = pos // 8
        bit_idx = pos % 8
        self.chromosomes[chrom][byte_idx] |= (1 << bit_idx)

    def mark_region(self, chrom: str, start: int, end: int):
        """
        Mark a region as covered (inclusive).

        Args:
            chrom: Chromosome name
            start: Start position (0-indexed)
            end: End position (inclusive)
        """
        if chrom not in self.chromosomes:
            return

        chrom_len = self.chrom_lengths[chrom]
        start = max(0, start)
        end = min(end, chrom_len - 1)

        if start > end:
            return

        # Mark each position
        # For efficiency, we can optimize this for large regions
        # but for now, simple iteration works
        bit_array = self.chromosomes[chrom]

        for pos in range(start, end + 1):
            byte_idx = pos // 8
            bit_idx = pos % 8
            bit_array[byte_idx] |= (1 << bit_idx)

    def mark_regions_batch(self, chrom: str, regions: List[Tuple[int, int]]):
        """
        Mark multiple regions as covered.

        Uses a fast approach: create boolean array, mark regions, pack to bits.

        Args:
            chrom: Chromosome name
            regions: List of (start, end) tuples (inclusive)
        """
        if chrom not in self.chromosomes:
            return

        if not regions:
            return

        chrom_len = self.chrom_lengths[chrom]
        bit_array = self.chromosomes[chrom]

        # Convert regions to numpy array for vectorized operations
        regions_arr = np.array(regions, dtype=np.int64)

        # Clip to valid range
        regions_arr[:, 0] = np.clip(regions_arr[:, 0], 0, chrom_len - 1)
        regions_arr[:, 1] = np.clip(regions_arr[:, 1], 0, chrom_len - 1)

        # Filter out invalid regions
        valid = regions_arr[:, 0] <= regions_arr[:, 1]
        regions_arr = regions_arr[valid]

        if len(regions_arr) == 0:
            return

        # Fast approach: use cumsum trick to mark regions
        # Create delta array: +1 at start, -1 at end+1
        # Then cumsum gives coverage depth, >0 means covered
        delta = np.zeros(chrom_len + 1, dtype=np.int32)
        np.add.at(delta, regions_arr[:, 0], 1)
        np.add.at(delta, regions_arr[:, 1] + 1, -1)

        # Cumsum and convert to boolean
        covered = np.cumsum(delta[:-1]) > 0

        # Pad to multiple of 8
        pad_len = (8 - (chrom_len % 8)) % 8
        if pad_len > 0:
            covered = np.concatenate([covered, np.zeros(pad_len, dtype=bool)])

        # Pack bits (numpy packbits uses MSB first, we want LSB first)
        # Reshape to (N, 8) and manually pack
        covered_reshaped = covered.reshape(-1, 8)
        powers = np.array([1, 2, 4, 8, 16, 32, 64, 128], dtype=np.uint8)
        new_bits = (covered_reshaped * powers).sum(axis=1).astype(np.uint8)

        # OR with existing bits
        bit_array[:len(new_bits)] |= new_bits

    def is_covered(self, chrom: str, pos: int) -> bool:
        """Check if a position is covered."""
        if chrom not in self.chromosomes:
            return False

        if pos < 0 or pos >= self.chrom_lengths[chrom]:
            return False

        byte_idx = pos // 8
        bit_idx = pos % 8
        return bool(self.chromosomes[chrom][byte_idx] & (1 << bit_idx))

    def count_covered(self, chrom: Optional[str] = None) -> int:
        """
        Count total covered positions.

        Args:
            chrom: If specified, count only for this chromosome.
                   If None, count for all chromosomes.

        Returns:
            Number of covered positions
        """
        if chrom is not None:
            if chrom not in self.chromosomes:
                return 0
            return self._count_bits(self.chromosomes[chrom])

        total = 0
        for chrom_array in self.chromosomes.values():
            total += self._count_bits(chrom_array)
        return total

    def _count_bits(self, arr: np.ndarray) -> int:
        """Count set bits in a numpy array using efficient popcount."""
        # Use numpy's unpackbits for efficient bit counting
        return int(np.unpackbits(arr).sum())

    def get_total_length(self) -> int:
        """Get total length of all chromosomes."""
        return sum(self.chrom_lengths.values())

    def get_coverage_stats(self) -> dict:
        """
        Get coverage statistics.

        Returns:
            Dictionary with coverage statistics
        """
        total_length = self.get_total_length()
        covered = self.count_covered()

        return {
            'total_length': total_length,
            'covered_bp': covered,
            'uncovered_bp': total_length - covered,
            'coverage_pct': (covered / total_length * 100) if total_length > 0 else 0,
            'chromosomes': {
                chrom: {
                    'length': self.chrom_lengths[chrom],
                    'covered': self.count_covered(chrom),
                }
                for chrom in self.chromosomes
            }
        }

    def count_new_coverage(self, chrom: str, regions: List[Tuple[int, int]]) -> int:
        """
        Count how many new positions would be covered by given regions.

        Args:
            chrom: Chromosome name
            regions: List of (start, end) tuples

        Returns:
            Number of positions not already covered
        """
        if chrom not in self.chromosomes:
            return 0

        if not regions:
            return 0

        chrom_len = self.chrom_lengths[chrom]
        bit_array = self.chromosomes[chrom]

        # Convert regions to numpy array
        regions_arr = np.array(regions, dtype=np.int64)

        # Clip to valid range
        regions_arr[:, 0] = np.clip(regions_arr[:, 0], 0, chrom_len - 1)
        regions_arr[:, 1] = np.clip(regions_arr[:, 1], 0, chrom_len - 1)

        # Filter out invalid regions
        valid = regions_arr[:, 0] <= regions_arr[:, 1]
        regions_arr = regions_arr[valid]

        if len(regions_arr) == 0:
            return 0

        # Use cumsum trick to find covered positions
        delta = np.zeros(chrom_len + 1, dtype=np.int32)
        np.add.at(delta, regions_arr[:, 0], 1)
        np.add.at(delta, regions_arr[:, 1] + 1, -1)

        # Get boolean array of new positions to cover
        new_covered = np.cumsum(delta[:-1]) > 0

        # Unpack existing coverage to boolean array
        existing_covered = np.unpackbits(bit_array).astype(bool)[:chrom_len]

        # Count positions that are new (would be covered but aren't already)
        new_only = new_covered & ~existing_covered

        return int(new_only.sum())

    def copy(self) -> 'CoverageTracker':
        """Create a deep copy of this tracker."""
        new_tracker = CoverageTracker()
        new_tracker.chrom_lengths = self.chrom_lengths.copy()
        new_tracker.chromosomes = {
            chrom: arr.copy() for chrom, arr in self.chromosomes.items()
        }
        return new_tracker
