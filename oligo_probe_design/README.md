# Oligo Probe Design

A Python package for designing specific oligonucleotide probes for single molecule RNA FISH. This package provides a local, self-contained implementation of the probe design algorithm originally developed in the Raj Lab.

## Features

- **Local execution**: No need for external servers or web services
- **Comprehensive masking**: Pseudogene, genome, GC run, and GC composition masking
- **Thermodynamic optimization**: RNA-DNA hybridization energetics using SantaLucia parameters
- **Flexible parameters**: Customizable oligo length, target melting temperature, and masking options
- **Multiple species support**: Human, mouse, fly, worm, rat, and cow
- **Command-line interface**: Easy to use from the command line
- **Python API**: Programmatic access for integration into workflows

## Installation

### Prerequisites

1. **Python 3.6+**
2. **Bowtie 1** (not Bowtie 2) - see installation instructions below

### Installing Bowtie 1

Bowtie 1 is required for sequence alignment and masking. **Note: This package requires Bowtie 1, not Bowtie 2.**

#### Option 1: Download pre-compiled binary (Recommended)

1. Download Bowtie 1 from: https://sourceforge.net/projects/bowtie-bio/files/bowtie/
2. Choose the latest version (e.g., bowtie-1.3.1)
3. Download the appropriate binary for your system
4. Extract and place in one of these locations:
   - `/usr/bin/bowtie` (system-wide, requires admin)
   - `$HOME/bowtie/bowtie` (recommended for most users)
   - `$HOME/Downloads/bowtie/bowtie`

#### Option 2: Set BOWTIEHOME environment variable

If you install Bowtie 1 in a custom location, set the environment variable:

```bash
export BOWTIEHOME=/path/to/your/bowtie-directory
```

Add this line to your `~/.bashrc` or `~/.bash_profile` to make it permanent.

#### Verify Bowtie 1 installation

```bash
# Test bowtie is found and working
bowtie --version
# Should show: bowtie version 1.x.x
```

### Installing Genome Indexes

You need to download pre-built bowtie indexes for your species of interest.

#### Option 1: Download from Raj Lab Dropbox (if available)

Contact the Raj Lab for access to pre-built indexes.

#### Option 2: Download from Bowtie website

1. Visit: http://bowtie-bio.sourceforge.net/tutorial.shtml#preb
2. Download indexes for your species (e.g., human, mouse)
3. Create an `indexes` directory in your bowtie installation:
   ```bash
   mkdir -p $HOME/bowtie/indexes
   ```
4. Extract indexes to this directory

#### Option 3: Create indexes directory in package

```bash
# Create indexes directory in the package
mkdir -p /path/to/ProbeDesign/oligo_probe_design/indexes/bowtie/indexes

# Download and extract genome indexes here
# Example for human:
cd /path/to/ProbeDesign/oligo_probe_design/indexes/bowtie/indexes
wget http://bowtie-bio.sourceforge.net/tutorial.shtml#preb
# Extract human genome, pseudogene, and mitochondrial indexes
```

### Required Indexes by Species

For each species, you need these bowtie indexes:

- **human**: `human`, `humanPseudo`, `humanMito`
- **mouse**: `mouse`, `mousePseudo`, `mouseMito`  
- **drosophila**: `drosophila`, `drosophilaPseudo`, `drosophilaMito`
- **elegans**: `celegans`, `celegansPseudo`, `celegansMito`
- **rat**: `rat`, `ratPseudo`
- **cow**: `cow`

## Usage

### Command Line Interface

#### Basic usage

```bash
# Design 48 probes for a human sequence
python -m oligo_probe_design input.fasta

# Design 32 probes for mouse
python -m oligo_probe_design input.fasta --n-oligos 32 --species mouse

# Custom thermodynamic parameters
python -m oligo_probe_design input.fasta --target-tm -25 --target-range -28 -22

# Disable specific masking options
python -m oligo_probe_design input.fasta --no-pseudogene-mask --no-gc-mask

# Specify output prefix
python -m oligo_probe_design input.fasta --output my_probes
```

#### Command line options

```
positional arguments:
  input_file            Input FASTA file containing target sequence

optional arguments:
  -h, --help            show this help message and exit
  --n-oligos N          Number of oligos to design (default: 48)
  --oligo-length N      Length of each oligo (default: 20)
  --spacer-length N     Minimum spacer between oligos (default: 2)
  --target-tm TM        Target free energy in kcal/mol (default: -23.0)
  --target-range MIN MAX
                        Allowable free energy range (default: -26.0 -20.0)
  --species {human,mouse,drosophila,elegans,rat,cow}
                        Species for masking databases (default: human)
  --no-pseudogene-mask  Disable pseudogene masking
  --no-genome-mask      Disable genome masking
  --no-gc-run-mask      Disable GC run masking
  --no-gc-mask          Disable GC composition masking
  --repeat-mask         Enable repeat masking (not yet implemented)
  --output PREFIX       Output file prefix
  --quiet               Suppress progress messages
```

### Python API

```python
from oligo_probe_design import design_probes

# Basic usage
results = design_probes('input.fasta')

# Custom parameters
results = design_probes(
    'input.fasta',
    n_oligos=32,
    species='mouse',
    target_tm=-25.0,
    target_range=[-28.0, -22.0],
    pseudogene_mask=True,
    genome_mask=True,
    output_file='my_probes'
)

# Access results
oligos = results['oligos']          # List of oligo sequences
positions = results['positions']    # Positions in input sequence
scores = results['scores']          # Quality scores
```

### Output Files

The program generates two output files:

1. **`{prefix}_oligos.txt`**: Tab-separated file with probe names and sequences
   ```
   Probe_01    TGATCGATCGATCGATCGAT
   Probe_02    CGATCGATCGATCGATCGAT
   ...
   ```

2. **`{prefix}_alignment.txt`**: Visual alignment showing probe positions on the target sequence

## Algorithm Details

The probe design algorithm implements the following steps:

1. **Sequence Loading**: Parse input FASTA file
2. **Masking**: Apply various masking strategies:
   - **Pseudogene masking**: Avoid regions similar to pseudogenes (16-mer, 0 mismatches)
   - **Genome masking**: Avoid highly repetitive genomic regions (multiple seed lengths)
   - **GC run masking**: Avoid runs of G or C nucleotides (≥7 in a row)
   - **GC composition masking**: Avoid oligos with extreme GC content
3. **Thermodynamic scoring**: Calculate RNA-DNA hybridization free energy for all possible oligos
4. **Optimization**: Use dynamic programming to find optimal oligo positions that minimize thermodynamic score while respecting spacing constraints

The algorithm matches the behavior of the original MATLAB `findprobesHD` function and produces equivalent results.

## Troubleshooting

### Common Issues

1. **"Could not find the bowtie executable!"**
   - Install Bowtie 1 (not Bowtie 2) following the instructions above
   - Check that `bowtie --version` works
   - Set `BOWTIEHOME` environment variable if needed

2. **"Provided 'species' database name 'X' is not known"**
   - Download the required bowtie indexes for your species
   - Check that index files are in the correct location
   - Verify index file names match expected pattern

3. **"No suitable oligos were found"**
   - Try relaxing thermodynamic constraints (`--target-range`)
   - Disable some masking options (`--no-pseudogene-mask`, etc.)
   - Check input sequence quality and length

4. **Alignment failures**
   - Ensure bowtie indexes are properly installed
   - Check network connectivity if downloading indexes
   - Verify index files are not corrupted

### Getting Help

For issues specific to this package:
1. Check that all dependencies are properly installed
2. Verify your input FASTA file is valid
3. Try running with `--quiet` disabled to see detailed progress

For issues with the original algorithm or biological interpretation:
- Refer to the original Raj Lab publications
- Check the MATLAB implementation in the parent directory

## Differences from MATLAB Version

This Python implementation:
- ✅ Implements the V4 algorithm with RNA-DNA thermodynamics
- ✅ Supports all major masking strategies
- ✅ Produces equivalent results to MATLAB version
- ❌ Does not yet implement repeat masking (RepeatMasker integration)
- ❌ Does not include the GUI interface
- ❌ Does not support all legacy parameter combinations

## Development

### Package Structure

```
oligo_probe_design/
├── __init__.py           # Package initialization
├── __main__.py           # CLI entry point  
├── cli.py                # Command line interface
├── README.md             # This file
└── core/                 # Core algorithms
    ├── __init__.py
    ├── fasta.py          # FASTA file handling
    ├── seq_utils.py      # Sequence utilities
    ├── bowtie_alignment.py # Bowtie integration
    ├── probe_design.py   # Core algorithms
    └── probe_designer.py # Main interface
```

## License

This package implements algorithms originally developed in the Raj Lab. Please cite the original publications when using this software for research.

## Citation

If you use this software, please cite:

[Original Raj Lab publications - to be added]