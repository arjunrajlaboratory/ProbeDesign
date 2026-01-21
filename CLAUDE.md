# Claude Code Context for ProbeDesign

This file provides context for Claude Code when working on this project.

## Project Overview

ProbeDesign is a tool for designing oligonucleotide probes for single molecule RNA FISH experiments. It was originally written in MATLAB and has been ported to Python for easier installation and use.

**Key goal**: The Python implementation should produce **identical results** to the MATLAB reference implementation (`probedesign/findprobesLocal.m`).

## Architecture

### Python Package (`src/probedesign/`)

| File | Purpose |
|------|---------|
| `cli.py` | Click-based command-line interface |
| `core.py` | Main probe design algorithm (badness calculation, DP optimization) |
| `thermodynamics.py` | Gibbs free energy and melting temperature calculations (Sugimoto 1995 params) |
| `masking.py` | Bowtie-based sequence masking (pseudogene, genome) |
| `sequence.py` | Sequence utilities (reverse complement, GC%, validation) |
| `fasta.py` | FASTA file I/O with junction marker handling |
| `output.py` | Output file generation (_oligos.txt, _seq.txt) |

### MATLAB Code (`probedesign/`)

| File | Purpose |
|------|---------|
| `findprobesLocal.m` | Main MATLAB function - **reference implementation** |
| `thermo_RNA_DNA.m` | Thermodynamics (Sugimoto 1995 parameters) |
| `find_best_matches.m` | Dynamic programming algorithm |
| `pseudogeneDBs/` | Pseudogene FASTA files for masking |

## Algorithm

1. **Read input**: Parse FASTA, concatenate multi-entry files with `>` junction markers
2. **Calculate badness**: For each position, compute `(gibbs - target)^2`. Positions with invalid chars or out-of-range Gibbs get `inf` badness.
3. **Apply masking**:
   - R mask: Repeat regions (from N's in input or separate repeatmask file)
   - P mask: Pseudogene alignments (bowtie to pseudogene DBs)
   - B mask: Genome alignments (bowtie to genome with multiple stringencies)
   - F mask: Thermodynamic filtering (badness == inf)
4. **Dynamic programming**: Find optimal probe placements minimizing average badness
5. **Generate output**: Create oligo and sequence alignment files

## Output Format

### `_oligos.txt`
Tab-separated: `index  GC%  Tm  Gibbs  sequence  name`

### `_seq.txt`
Visual alignment showing:
- Line 1: Original sequence
- Line 2: R mask (if repeat masking enabled)
- Line 3+: P mask, B mask (if enabled)
- Last mask line: F mask (thermodynamic filtering)
- Probe annotations below each block

## Test Cases (`test_cases/`)

| Test Case | Description | Expected Match |
|-----------|-------------|----------------|
| `KRT19_withUTRs/` | 6 probes, pseudogene masking | 100% |
| `CDKN1A_32/` | 32 probes, repeat masking | 100% (with genomemaskoff) |
| `EIF1_CDS_HCR/` | HCR probes (longer oligos) | Validated |

### Running Test Validation

```bash
# Basic test (no masking)
probedesign design test_cases/KRT19_withUTRs/KRT19_withUTRs.fa --probes 6

# With repeat masking (manual)
probedesign design test_cases/CDKN1A_32/CDKN1A.fa --probes 32 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa

# With bowtie masking (requires bowtie + indexes)
BOWTIEHOME=/opt/homebrew/Caskroom/miniforge/base/bin probedesign design \
  test_cases/KRT19_withUTRs/KRT19_withUTRs.fa --probes 6 --pseudogene-mask
```

## Important Implementation Details

### Thermodynamics

- Uses Sugimoto 1995 RNA-DNA hybrid parameters
- Salt concentration: 0.33 M
- Primer concentration: 50 µM
- Default target Gibbs: -23 kcal/mol (range: -26 to -20)

### F Mask Logic

The F mask shows thermodynamic filtering - positions where a probe **cannot start**:
- `F` at position i if `badness[i] == inf` OR `i >= goodlen` (past valid positions)
- Sequence char at position i if `badness[i]` is finite (valid probe start)

This matches MATLAB `mask_string(inseq, badness==inf, 'F')`.

### Junction Handling

Multi-entry FASTA files use `>` to mark exon junctions. Probes cannot span junctions.

### Repeat Masking

Two modes:
1. **Auto-detect**: N's in input file are treated as masked regions
2. **Manual file**: `--repeatmask-file` provides separate file with N's marking repeats

## Bowtie Setup

See [BOWTIE.md](BOWTIE.md) for installation instructions.

Key points:
- Uses Bowtie 1 (not Bowtie 2) for short read alignment
- Genome indexes: `bowtie_indexes/` directory
- Pseudogene DBs: `probedesign/pseudogeneDBs/`
- Set `BOWTIEHOME` env var if bowtie not in PATH

## Common Development Tasks

### Adding a new CLI option

1. Add Click option in `cli.py`
2. Pass parameter to `design_probes()` in `core.py`
3. Update `design_probes()` function signature and logic
4. Update README.md and this file

### Validating against MATLAB

1. Run MATLAB command with specific options
2. Save `_oligos.txt` and `_seq.txt` output
3. Run Python with equivalent options
4. Compare outputs with `diff`

### Debugging probe differences

1. Check F mask line matches (thermodynamic filtering)
2. Check mask strings (R, P, B) match
3. Verify badness calculation for specific positions
4. Compare probe positions and sequences

## Code Style

- Python 3.8+ with type hints
- Click for CLI
- No external dependencies except numpy (optional)
- Match MATLAB output format exactly for validation
