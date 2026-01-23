# ProbeDesign

Single molecule RNA FISH probe design software from the Raj Lab.

## Overview

ProbeDesign creates oligonucleotide probes for single molecule RNA FISH experiments. It selects optimal probe sequences based on:

- **Thermodynamic properties**: Targets specific Gibbs free energy values for RNA-DNA hybridization
- **Sequence masking**: Avoids cross-hybridization by masking regions that align to pseudogenes or repetitive genomic sequences
- **Optimal spacing**: Uses dynamic programming to find the best probe placements

## Installation

### Python CLI (Recommended)

```bash
# Clone the repository
git clone https://github.com/arjunrajlab/ProbeDesign.git
cd ProbeDesign

# Install with pip (requires Python 3.8+)
pip install -e .

# Verify installation
probedesign --version
```

### Bowtie Setup (for masking)

Masking requires [Bowtie 1](http://bowtie-bio.sourceforge.net/) and genome indexes. See [BOWTIE.md](BOWTIE.md) for detailed setup instructions.

**Important**: Install bowtie via conda/mamba (NOT homebrew, which has a different tool with the same name).

```bash
# Install bowtie via conda/mamba
mamba install bowtie

# Verify it's the bioinformatics bowtie
bowtie --version  # Should show "bowtie-align-s version 1.x.x"

# Create indexes directory
mkdir -p bowtie_indexes && cd bowtie_indexes

# Build pseudogene index from included FASTA files
bowtie-build ../probedesign/pseudogeneDBs/human.fasta humanPseudo

# Download human genome index (~3GB)
curl -O ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/GRCh38_no_alt.zip
unzip GRCh38_no_alt.zip
```

For other species, build the corresponding pseudogene index:
```bash
bowtie-build ../probedesign/pseudogeneDBs/mouse.fasta mousePseudo
bowtie-build ../probedesign/pseudogeneDBs/elegans.fasta celegansPseudo
bowtie-build ../probedesign/pseudogeneDBs/drosophila.fasta drosophilaPseudo
```

## Quick Start

```bash
# Basic usage (no masking)
probedesign design input.fa

# With pseudogene and genome masking (recommended)
probedesign design input.fa --pseudogene-mask --genome-mask --index-dir bowtie_indexes

# With automatic repeat masking (requires RepeatMasker installed locally)
probedesign design input.fa --repeatmask

# With manual repeat masking (using RepeatMasker output file)
probedesign design input.fa --repeatmask-file input_repeatmasked.fa

# Full example with all options
probedesign design input.fa \
  --probes 48 \
  --oligo-length 20 \
  --spacer-length 2 \
  --target-gibbs -23 \
  --allowable-gibbs -26,-20 \
  --repeatmask \
  --pseudogene-mask \
  --genome-mask \
  --species human \
  --index-dir bowtie_indexes \
  --output MyProbes
```

## CLI Reference

```
probedesign design [OPTIONS] INPUT_FILE

Arguments:
  INPUT_FILE                FASTA file containing target sequence

Options:
  -n, --probes INTEGER      Number of probes to design (default: 48)
  -l, --oligo-length INT    Length of each oligonucleotide (default: 20)
  -s, --spacer-length INT   Minimum gap between probes (default: 2)
  -g, --target-gibbs FLOAT  Target Gibbs free energy in kcal/mol (default: -23)
  --allowable-gibbs TEXT    Allowable Gibbs FE range as min,max (default: -26,-20)
  -o, --output TEXT         Output file prefix (default: input filename)
  -q, --quiet               Suppress output to stdout
  --species [human|mouse|elegans|drosophila|rat]
                            Species for masking databases (default: human)
  --repeatmask              Run RepeatMasker automatically (requires local install)
  --repeatmask-file PATH    FASTA file with N's marking repeat regions (manual)
  --pseudogene-mask         Mask regions that align to pseudogenes
  --genome-mask             Mask repetitive genomic regions
  --index-dir PATH          Directory containing bowtie indexes
```

## Output Files

ProbeDesign generates two output files:

### `<output>_oligos.txt`

Tab-separated file with probe details:

```
Index  GC%  Tm    Gibbs   Sequence              Name
1      55   66.2  -24.7   gtctcagaagctgcgattcg  MyProbes_1
2      65   63.2  -26.0   gaatgctgggcgcgcgaaag  MyProbes_2
...
```

### `<output>_seq.txt`

Visual alignment of probes against the input sequence, showing masked regions and probe positions.

## Masking Options

### Pseudogene Masking (`--pseudogene-mask`)

Masks regions that align to known pseudogene sequences. Uses 16-mer alignment with bowtie. Short runs (<18bp) are removed to avoid spurious masking.

### Genome Masking (`--genome-mask`)

Masks repetitive genomic regions using multiple alignment stringencies:
- 12-mer with threshold 4000 hits (very repetitive)
- 14-mer with threshold 500 hits
- 16-mer with threshold 20 hits

### Automatic Repeat Masking (`--repeatmask`)

If you have RepeatMasker installed locally, ProbeDesign can run it automatically:

```bash
probedesign design input.fa --repeatmask --probes 32
```

**Quick Setup (Human/Mouse):**
```bash
# Install RepeatMasker
mamba install -c bioconda -c conda-forge repeatmasker

# Download Mammalia database (~8.9GB compressed, ~56GB extracted)
cd /opt/homebrew/Caskroom/miniforge/base/share/RepeatMasker/Libraries/famdb
curl -O https://www.dfam.org/releases/current/families/FamDB/dfam39_full.7.h5.gz
gunzip dfam39_full.7.h5.gz
```

**Note:** The database requires ~65GB disk space. See [REPEATMASKER.md](REPEATMASKER.md) for detailed setup instructions, troubleshooting, and the web-based alternative if disk space is limited.

### Manual Repeat Masking (`--repeatmask-file`)

For sequences with repetitive elements, you can provide a pre-masked file:

**Option A**: Use RepeatMasker's website
1. Upload your FASTA to [RepeatMasker's website](http://www.repeatmasker.org/)
2. Download the masked output file (contains N's where repeats are)
3. Use the `--repeatmask-file` option:

**Option B**: Run RepeatMasker locally
1. Run: `RepeatMasker -species human input.fa`
2. Use the output: `probedesign design input.fa --repeatmask-file input.fa.masked`

```bash
probedesign design input.fa --repeatmask-file input_repeatmasked.fa --output MyProbes
```

**Note**: You cannot use both `--repeatmask` and `--repeatmask-file` together.

### Species Support

Masking databases are available for: `human`, `mouse`, `elegans`, `drosophila`, `rat`

## Thermodynamics

ProbeDesign uses Sugimoto 1995 parameters for RNA-DNA hybrid thermodynamics:

- **Gibbs Free Energy**: Calculated at 37°C
- **Melting Temperature**: Calculated with salt=0.33M, primer concentration=50µM
- **Default target**: -23 kcal/mol (range: -26 to -20)

For HCR probes (longer oligos), use different parameters:
```bash
probedesign design input.fa -l 52 --target-gibbs -60 --allowable-gibbs -80,-40
```

## MATLAB Version (Legacy)

The original MATLAB implementation is in `probedesign/findprobesLocal.m`. It requires MATLAB and has additional features like repeat masking via RepeatMasker.

```matlab
% Basic usage
findprobesLocal('input.fa', 48)

% With options
findprobesLocal('input.fa', 32, ...
    'outfilename', 'MyProbes', ...
    'pseudogenemask', true, ...
    'genomemask', true, ...
    'repeatmask', false, ...
    'species', 'human')
```

### RepeatMasker Integration

**Python CLI (recommended):**

Option 1 - Automatic (requires local RepeatMasker):
```bash
probedesign design mCherry.fasta --repeatmask --probes 32
```

Option 2 - Manual (using pre-masked file):
```bash
probedesign design mCherry.fasta --repeatmask-file mCherry_repeatmasked.fa --output mCherry_probes
```

See [REPEATMASKER.md](REPEATMASKER.md) for RepeatMasker installation instructions.

**MATLAB:**
```matlab
findprobesLocal('mCherry.fasta', 32, ...
    'repeatmask', true, ...
    'repeatmaskmanual', true, ...
    'repeatmaskfile', 'mCherry_repeatmasked.fa', ...
    'species', 'mouse')
```

## Testing

See [TEST.md](TEST.md) for validation instructions and test cases.

## Project Structure

```
ProbeDesign/
├── src/probedesign/          # Python package
│   ├── cli.py                # Command-line interface
│   ├── core.py               # Main probe design algorithm
│   ├── thermodynamics.py     # Gibbs FE and Tm calculations
│   ├── masking.py            # Bowtie-based sequence masking
│   ├── sequence.py           # Sequence utilities
│   ├── fasta.py              # FASTA I/O
│   └── output.py             # Output file generation
├── probedesign/              # MATLAB code
│   ├── findprobesLocal.m     # Main MATLAB function
│   └── pseudogeneDBs/        # Pseudogene FASTA files
├── DesignServer/             # Legacy Python server code
├── test_cases/               # Validation test cases
├── bowtie_indexes/           # Genome indexes (gitignored)
├── pyproject.toml            # Python package config
├── BOWTIE.md                 # Bowtie setup instructions
├── REPEATMASKER.md           # RepeatMasker setup instructions
└── TEST.md                   # Testing documentation
```

## Algorithm

1. **Read input**: Parse FASTA file, clean sequence
2. **Calculate badness**: For each position, compute squared distance from target Gibbs FE
3. **Apply masking**: Mark regions aligning to pseudogenes/repetitive sequences as infinite badness
4. **Dynamic programming**: Find optimal probe placements minimizing average badness
5. **Generate output**: Create oligo and sequence alignment files

## License

MIT License - See LICENSE file for details.

## Citation

If you use this software, please cite the Raj Lab.
