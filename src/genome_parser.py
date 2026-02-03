"""Stream-parse FASTA files for genome sequences."""

import gzip
from pathlib import Path
from typing import Generator, Tuple, Set

# Primary chromosomes to process
PRIMARY_CHROMOSOMES = {
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
    'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
    'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
    'chrX', 'chrY', 'chrM',
    # Also support without 'chr' prefix
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
    '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
    'X', 'Y', 'M', 'MT',
}


def normalize_chrom_name(name: str) -> str:
    """Normalize chromosome name to standard format (chr1, chr2, etc.)."""
    name = name.strip()
    # Handle names like ">chr1 extra info"
    name = name.split()[0]
    if name.startswith('>'):
        name = name[1:]

    # Add 'chr' prefix if missing
    if not name.startswith('chr'):
        if name == 'MT':
            return 'chrM'
        return f'chr{name}'
    return name


def is_primary_chromosome(name: str) -> bool:
    """Check if chromosome name is a primary chromosome."""
    # Extract just the chromosome identifier
    name = name.strip()
    if name.startswith('>'):
        name = name[1:]
    name = name.split()[0]

    return name in PRIMARY_CHROMOSOMES


def parse_fasta(filepath: str) -> Generator[Tuple[str, str], None, None]:
    """
    Parse a FASTA file and yield (chromosome_name, sequence) tuples.

    Supports gzipped files (detected by .gz extension).
    Only yields primary chromosomes (chr1-22, chrX, chrY, chrM).

    Args:
        filepath: Path to FASTA file (can be .gz compressed)

    Yields:
        Tuples of (chromosome_name, sequence)
    """
    filepath = Path(filepath)

    # Open file (handle gzip)
    if filepath.suffix == '.gz':
        file_handle = gzip.open(filepath, 'rt')
    else:
        file_handle = open(filepath, 'r')

    try:
        current_chrom = None
        current_seq = []
        include_current = False

        for line in file_handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                # Yield previous chromosome if we have one
                if current_chrom is not None and include_current:
                    yield (current_chrom, ''.join(current_seq).upper())

                # Start new chromosome
                if is_primary_chromosome(line):
                    current_chrom = normalize_chrom_name(line)
                    current_seq = []
                    include_current = True
                else:
                    include_current = False
                    current_chrom = None
                    current_seq = []
            else:
                if include_current:
                    current_seq.append(line)

        # Yield last chromosome
        if current_chrom is not None and include_current:
            yield (current_chrom, ''.join(current_seq).upper())

    finally:
        file_handle.close()


def get_genome_stats(filepath: str) -> dict:
    """
    Get statistics about the genome file.

    Args:
        filepath: Path to FASTA file

    Returns:
        Dictionary with genome statistics
    """
    total_length = 0
    chromosomes = {}

    for chrom, seq in parse_fasta(filepath):
        chrom_len = len(seq)
        chromosomes[chrom] = chrom_len
        total_length += chrom_len

    return {
        'total_length': total_length,
        'chromosomes': chromosomes,
        'num_chromosomes': len(chromosomes),
    }
