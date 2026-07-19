# BetterCrispr

A command-line tool that measures what fraction of a genome is targetable by CRISPR for a given set of PAM sequences.

## What it does

Given a genome FASTA file and one or more PAM patterns (e.g. `NGG` for SpCas9), BetterCrispr scans both strands, locates every PAM site, projects the editable region around each Cas9 cut site, and reports the total fraction of the genome that can be reached. It can also **optimize**: given the PAMs you already have, it searches all PAMs of a chosen length and ranks the ones that would add the most *new* coverage.

Three subcommands are available:

- `coverage` — report the combined coverage for one or more PAMs.
- `optimize` — rank the PAMs that add the most additional coverage on top of your current set.
- `analyze` — an interactive REPL to add PAMs, inspect stats, optimize, and export results incrementally.

PAM patterns accept [IUPAC ambiguity codes](https://en.wikipedia.org/wiki/Nucleic_acid_notation) (`N`, `R`, `Y`, `W`, `V`, etc.), so `NGG`, `TTTN`, and `NRG` are all valid.

### Coverage models

The `--coverage-mode` flag controls how many base pairs each cut site is counted as covering:

| Mode        | Meaning                                                        |
|-------------|---------------------------------------------------------------|
| `cut_site`  | Only the exact cut position (1 bp).                           |
| `editing`   | ~21 bp HDR editing window around the cut site. **(default)**  |
| `guide`     | The full 20 bp guide-binding region (overstates real editing reach). |
| `prime`     | ~30 bp downstream of the nick (prime-editing window).         |

## Requirements

- **Python 3.8+**
- Python packages: `numpy` and `tqdm` (see `requirements.txt`).
- **A genome FASTA file that you supply yourself.** No genome ships with this repo — the analysis is meaningless without one. Any FASTA works; a plain or gzipped human reference such as **hg38** is the typical input.
  - Download example (UCSC hg38, ~950 MB compressed):
    ```bash
    mkdir -p genome
    curl -o genome/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    ```
  - Only **primary chromosomes** are analyzed: `chr1`–`chr22`, `chrX`, `chrY`, `chrM` (with or without the `chr` prefix). Scaffolds, alt contigs, and unplaced sequences in the FASTA are ignored.
  - `.gz` files are read directly — no need to decompress first. FASTA files are git-ignored by design, so they never get committed.

## Install

```bash
git clone https://github.com/Mkrolick/BetterCrispr.git
cd BetterCrispr
python3 -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Usage

Run everything through `main.py`. The `--genome` flag is required by every subcommand.

### Analyze coverage for a PAM

```bash
python main.py coverage --pam NGG --genome genome/hg38.fa.gz
```

Expected output (numbers depend on your genome):

```
PAM Coverage Analysis
==================================================
Genome: genome/hg38.fa.gz
Genome length: 3,088,269,832 bp
Coverage mode: editing (Efficient HDR editing window (~21bp around cut site))

PAM Sites:
  NGG: 246,153,024 sites

Combined Coverage: 2,801,344,112 bp (90.71%)
Uncovered: 286,925,720 bp (9.29%)
```

Multiple PAMs are comma-separated and combined into one coverage figure:

```bash
python main.py coverage --pam NGG,TTTN --genome genome/hg38.fa.gz
```

Add `--json` for machine-readable output, and `--quiet` / `-q` to suppress the load-progress lines.

### Find the best next PAM

Given the PAMs you already use, rank the PAMs of a chosen length that add the most new coverage:

```bash
python main.py optimize --current-pams NGG --pam-length 3 --top-n 10 --genome genome/hg38.fa.gz
```

```
Top 10 PAMs to Maximize Coverage:
--------------------------------------------------
   1. TTTN: +... bp (+...%) → ...% total
   ...
```

`--current-pams` is optional (omit it to rank from a blank slate). `--pam-length` defaults to `3` and `--top-n` to `10`.

### Interactive mode

```bash
python main.py analyze --genome genome/hg38.fa.gz
```

Then type commands at the `>>>` prompt:

```
add NGG          # add a PAM and update coverage
stats            # show current coverage statistics
optimize 3       # find the best 3-length PAMs to add next (optionally: optimize 3 5)
export out.json  # write current results to a JSON file
quit             # exit
```

### Common flags

| Flag              | Applies to           | Description                                                  |
|-------------------|----------------------|-------------------------------------------------------------|
| `--genome PATH`   | all                  | Path to the genome FASTA (`.fa`, `.fasta`, `.fna`, or `.gz`). **Required.** |
| `--pam SEQ[,SEQ]` | `coverage`           | Comma-separated PAM pattern(s). **Required.**               |
| `--current-pams`  | `optimize`           | Comma-separated PAMs you already have (optional).           |
| `--pam-length N`  | `optimize`           | Length of candidate PAMs to evaluate (default 3).           |
| `--top-n N`       | `optimize`           | Number of ranked PAMs to return (default 10).               |
| `--coverage-mode` | `coverage`, `optimize`, `analyze` | `cut_site` \| `editing` \| `guide` \| `prime` (default `editing`). |
| `--json`          | `coverage`, `optimize` | Emit JSON instead of the text report.                     |
| `--quiet`, `-q`   | `coverage`, `optimize` | Suppress genome-loading progress output.                  |

## Notes

- **The whole genome is loaded into memory.** A full human genome needs several GB of RAM and a minute or two to parse; `optimize` with a large `--pam-length` multiplies the work because it evaluates every PAM of that length. Start with a small chromosome or a `--pam-length` of 2–3 to gauge runtime.
- The `coverage --json` output reports combined coverage only; the per-PAM `coverage_bp` field is intentionally `null` because individual PAM contributions are not tracked separately (the coverage bit-vector is shared).
- Cut-site geometry assumes SpCas9 (cut 3 bp upstream of the PAM). The `prime` mode is a coarse downstream-window approximation, not a validated prime-editing efficiency model.
