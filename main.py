#!/usr/bin/env python3
"""
CRISPR PAM Coverage Analyzer

Analyze what fraction of the human genome is targetable by CRISPR
given specific PAM sequences.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional

from tqdm import tqdm

from src.genome_parser import parse_fasta, get_genome_stats
from src.coverage import CoverageTracker
from src.optimizer import analyze_pam_coverage, analyze_multiple_pams, find_best_next_pam
from src.pam_finder import COVERAGE_MODES, COVERAGE_CUT_SITE, COVERAGE_EDITING_WINDOW, COVERAGE_GUIDE_BINDING


def load_genome(filepath: str, show_progress: bool = True) -> dict:
    """Load genome sequences into memory."""
    sequences = {}

    if show_progress:
        print(f"Loading genome from {filepath}...")

    for chrom, seq in parse_fasta(filepath):
        sequences[chrom] = seq
        if show_progress:
            print(f"  Loaded {chrom}: {len(seq):,} bp")

    if show_progress:
        total = sum(len(s) for s in sequences.values())
        print(f"Total: {total:,} bp across {len(sequences)} chromosomes\n")

    return sequences


def format_number(n: int) -> str:
    """Format large number with commas."""
    return f"{n:,}"


def cmd_coverage(args):
    """Handle the coverage command."""
    # Load genome
    sequences = load_genome(args.genome, show_progress=not args.quiet)

    # Parse PAMs
    pams = [p.strip().upper() for p in args.pam.split(',')]

    # Get coverage mode
    coverage_mode = args.coverage_mode

    # Analyze coverage
    pam_sites, coverage = analyze_multiple_pams(
        sequences, pams, coverage_mode=coverage_mode, show_progress=not args.quiet
    )

    # Get stats
    stats = coverage.get_coverage_stats()

    # Output results
    mode_desc = COVERAGE_MODES[coverage_mode]['description']

    if args.json:
        output = {
            'genome': str(args.genome),
            'genome_length': stats['total_length'],
            'coverage_mode': coverage_mode,
            'coverage_mode_description': mode_desc,
            'pams': {
                pam: {
                    'sites': pam_sites[pam],
                    'coverage_bp': None  # Individual PAM coverage not tracked separately
                }
                for pam in pams
            },
            'total_coverage_bp': stats['covered_bp'],
            'total_coverage_pct': round(stats['coverage_pct'], 4),
            'uncovered_bp': stats['uncovered_bp'],
        }
        print(json.dumps(output, indent=2))
    else:
        print("PAM Coverage Analysis")
        print("=" * 50)
        print(f"Genome: {args.genome}")
        print(f"Genome length: {format_number(stats['total_length'])} bp")
        print(f"Coverage mode: {coverage_mode} ({mode_desc})")
        print()
        print("PAM Sites:")
        for pam in pams:
            print(f"  {pam}: {format_number(pam_sites[pam])} sites")
        print()
        print(f"Combined Coverage: {format_number(stats['covered_bp'])} bp ({stats['coverage_pct']:.2f}%)")
        print(f"Uncovered: {format_number(stats['uncovered_bp'])} bp ({100 - stats['coverage_pct']:.2f}%)")


def cmd_optimize(args):
    """Handle the optimize command."""
    # Load genome
    sequences = load_genome(args.genome, show_progress=not args.quiet)

    # Get coverage mode
    coverage_mode = args.coverage_mode

    # Initialize coverage tracker
    coverage = CoverageTracker()
    for chrom, seq in sequences.items():
        coverage.add_chromosome(chrom, len(seq))

    pam_sites = {}

    # Analyze current PAMs if provided
    if args.current_pams:
        current_pams = [p.strip().upper() for p in args.current_pams.split(',')]
        pam_sites, coverage = analyze_multiple_pams(
            sequences, current_pams, coverage_mode=coverage_mode, show_progress=not args.quiet
        )

        if not args.quiet:
            stats = coverage.get_coverage_stats()
            print(f"Current coverage: {format_number(stats['covered_bp'])} bp ({stats['coverage_pct']:.2f}%)")
            print()

    # Find best next PAMs
    if not args.quiet:
        print(f"Finding best PAMs of length {args.pam_length}...")

    results = find_best_next_pam(
        sequences, coverage, args.pam_length, args.top_n,
        coverage_mode=coverage_mode, show_progress=not args.quiet
    )

    # Get current stats for calculating totals
    current_stats = coverage.get_coverage_stats()
    current_coverage_pct = current_stats['coverage_pct']
    total_length = current_stats['total_length']

    # Output results
    mode_desc = COVERAGE_MODES[coverage_mode]['description']

    if args.json:
        output = {
            'genome': str(args.genome),
            'genome_length': total_length,
            'coverage_mode': coverage_mode,
            'coverage_mode_description': mode_desc,
            'current_pams': args.current_pams.split(',') if args.current_pams else [],
            'current_coverage_bp': current_stats['covered_bp'],
            'current_coverage_pct': round(current_coverage_pct, 4),
            'recommendations': [
                {
                    'pam': pam,
                    'additional_bp': additional_bp,
                    'additional_pct': round(additional_pct, 4),
                    'total_coverage_pct': round(current_coverage_pct + additional_pct, 4),
                }
                for pam, additional_bp, additional_pct in results
            ]
        }
        print(json.dumps(output, indent=2))
    else:
        print()
        print(f"Top {args.top_n} PAMs to Maximize Coverage:")
        print("-" * 50)
        for i, (pam, additional_bp, additional_pct) in enumerate(results, 1):
            total_pct = current_coverage_pct + additional_pct
            print(f"  {i:2}. {pam}: +{format_number(additional_bp)} bp (+{additional_pct:.2f}%) → {total_pct:.2f}% total")


def cmd_analyze(args):
    """Handle the interactive analyze command."""
    # Load genome
    sequences = load_genome(args.genome, show_progress=True)

    # Get coverage mode
    coverage_mode = args.coverage_mode
    mode_desc = COVERAGE_MODES[coverage_mode]['description']

    # Initialize coverage
    coverage = CoverageTracker()
    for chrom, seq in sequences.items():
        coverage.add_chromosome(chrom, len(seq))

    all_pam_sites = {}
    stats = coverage.get_coverage_stats()

    print("\nInteractive PAM Analysis")
    print("=" * 50)
    print(f"Genome length: {format_number(stats['total_length'])} bp")
    print(f"Coverage mode: {coverage_mode} ({mode_desc})")
    print("\nCommands:")
    print("  add <PAM>     - Add a PAM (e.g., 'add NGG')")
    print("  stats         - Show current coverage statistics")
    print("  optimize <N>  - Find best N-length PAMs to add next")
    print("  export <file> - Export results to JSON file")
    print("  quit          - Exit")
    print()

    while True:
        try:
            cmd = input(">>> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nGoodbye!")
            break

        if not cmd:
            continue

        parts = cmd.split()
        action = parts[0].lower()

        if action == 'quit' or action == 'exit':
            print("Goodbye!")
            break

        elif action == 'add' and len(parts) >= 2:
            pam = parts[1].upper()
            try:
                sites, coverage = analyze_pam_coverage(
                    sequences, pam, coverage,
                    coverage_mode=coverage_mode, show_progress=True
                )
                all_pam_sites[pam] = sites
                stats = coverage.get_coverage_stats()
                print(f"Added {pam}: {format_number(sites)} sites")
                print(f"Coverage now: {format_number(stats['covered_bp'])} bp ({stats['coverage_pct']:.2f}%)")
            except Exception as e:
                print(f"Error: {e}")

        elif action == 'stats':
            stats = coverage.get_coverage_stats()
            print(f"\nCurrent Coverage Statistics")
            print("-" * 30)
            if all_pam_sites:
                print("PAMs analyzed:")
                for pam, sites in all_pam_sites.items():
                    print(f"  {pam}: {format_number(sites)} sites")
                print()
            print(f"Total coverage: {format_number(stats['covered_bp'])} bp ({stats['coverage_pct']:.2f}%)")
            print(f"Uncovered: {format_number(stats['uncovered_bp'])} bp ({100 - stats['coverage_pct']:.2f}%)")

        elif action == 'optimize' and len(parts) >= 2:
            try:
                pam_length = int(parts[1])
                top_n = int(parts[2]) if len(parts) > 2 else 10
                results = find_best_next_pam(
                    sequences, coverage, pam_length, top_n,
                    coverage_mode=coverage_mode, show_progress=True
                )
                current_pct = coverage.get_coverage_stats()['coverage_pct']
                print(f"\nTop {top_n} PAMs of length {pam_length}:")
                for i, (pam, additional_bp, additional_pct) in enumerate(results, 1):
                    total_pct = current_pct + additional_pct
                    print(f"  {i:2}. {pam}: +{format_number(additional_bp)} bp (+{additional_pct:.2f}%) → {total_pct:.2f}%")
            except ValueError:
                print("Usage: optimize <pam_length> [top_n]")

        elif action == 'export' and len(parts) >= 2:
            filepath = parts[1]
            stats = coverage.get_coverage_stats()
            output = {
                'genome': str(args.genome),
                'genome_length': stats['total_length'],
                'coverage_mode': coverage_mode,
                'coverage_mode_description': mode_desc,
                'pams': {
                    pam: {'sites': sites}
                    for pam, sites in all_pam_sites.items()
                },
                'total_coverage_bp': stats['covered_bp'],
                'total_coverage_pct': round(stats['coverage_pct'], 4),
            }
            with open(filepath, 'w') as f:
                json.dump(output, f, indent=2)
            print(f"Exported to {filepath}")

        else:
            print("Unknown command. Type 'quit' to exit.")


def main():
    parser = argparse.ArgumentParser(
        description='CRISPR PAM Coverage Analyzer',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze coverage for a PAM
  python main.py coverage --pam NGG --genome genome/hg38.fa.gz

  # Analyze multiple PAMs
  python main.py coverage --pam NGG,TTTN --genome genome/hg38.fa.gz

  # Find optimal next PAM given current PAMs
  python main.py optimize --current-pams NGG --pam-length 3 --genome genome/hg38.fa.gz

  # Interactive analysis
  python main.py analyze --genome genome/hg38.fa.gz
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Coverage command
    coverage_parser = subparsers.add_parser('coverage', help='Analyze PAM coverage')
    coverage_parser.add_argument('--pam', required=True,
                                 help='PAM sequence(s), comma-separated (e.g., NGG or NGG,TTTN)')
    coverage_parser.add_argument('--genome', required=True,
                                 help='Path to genome FASTA file (can be gzipped)')
    coverage_parser.add_argument('--coverage-mode', choices=['cut_site', 'editing', 'guide', 'prime'],
                                 default='editing',
                                 help='Coverage model: cut_site (1bp), editing (~21bp HDR window), '
                                      'guide (20bp binding region), prime (~30bp downstream). Default: editing')
    coverage_parser.add_argument('--json', action='store_true',
                                 help='Output in JSON format')
    coverage_parser.add_argument('--quiet', '-q', action='store_true',
                                 help='Suppress progress output')

    # Optimize command
    optimize_parser = subparsers.add_parser('optimize', help='Find optimal next PAM')
    optimize_parser.add_argument('--current-pams',
                                 help='Current PAM sequence(s), comma-separated')
    optimize_parser.add_argument('--pam-length', type=int, default=3,
                                 help='Length of PAM sequences to evaluate (default: 3)')
    optimize_parser.add_argument('--top-n', type=int, default=10,
                                 help='Number of top PAMs to return (default: 10)')
    optimize_parser.add_argument('--genome', required=True,
                                 help='Path to genome FASTA file')
    optimize_parser.add_argument('--coverage-mode', choices=['cut_site', 'editing', 'guide', 'prime'],
                                 default='editing',
                                 help='Coverage model: cut_site (1bp), editing (~21bp HDR window), '
                                      'guide (20bp binding region), prime (~30bp downstream). Default: editing')
    optimize_parser.add_argument('--json', action='store_true',
                                 help='Output in JSON format')
    optimize_parser.add_argument('--quiet', '-q', action='store_true',
                                 help='Suppress progress output')

    # Analyze command (interactive)
    analyze_parser = subparsers.add_parser('analyze', help='Interactive analysis mode')
    analyze_parser.add_argument('--genome', required=True,
                                help='Path to genome FASTA file')
    analyze_parser.add_argument('--coverage-mode', choices=['cut_site', 'editing', 'guide', 'prime'],
                                default='editing',
                                help='Coverage model: cut_site (1bp), editing (~21bp HDR window), '
                                     'guide (20bp binding region). Default: editing')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == 'coverage':
        cmd_coverage(args)
    elif args.command == 'optimize':
        cmd_optimize(args)
    elif args.command == 'analyze':
        cmd_analyze(args)


if __name__ == '__main__':
    main()
