# Bowtie Installation and Setup

This document describes how to install Bowtie and set up genome indexes for probe design masking (pseudogene mask, genome mask).

## Overview

Bowtie is an ultrafast, memory-efficient short read aligner used by ProbeDesign to identify sequences that align to pseudogenes or multiple genomic locations. These alignments are used to mask problematic regions during probe design.

**Important Notes**:
- This project uses **Bowtie 1** (not Bowtie 2). Bowtie 1 is optimized for short reads (≤50bp) which matches our oligo lengths.
- Install bowtie via **conda/mamba** (bioinformatics tool), NOT via homebrew (which has a different "bowtie" package)
- You can verify you have the right bowtie by running: `bowtie --version` - it should show "bowtie-align-s version 1.x.x"

## Installation

### Option 1: Using Conda/Mamba (Recommended)

1. **Install Miniforge** (if not already installed):
   ```bash
   brew install --cask miniforge
   conda init "$(basename "${SHELL}")"
   # Restart your terminal
   ```

2. **Configure Bioconda channels**:
   ```bash
   conda config --add channels defaults
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda config --set channel_priority strict
   ```

3. **Install Bowtie**:
   ```bash
   mamba install bowtie
   # or: conda install bowtie
   ```

4. **Verify installation**:
   ```bash
   conda activate base
   bowtie --version
   # Should show: bowtie-align-s version 1.3.1 (or similar)
   ```

### Option 2: Direct Download

Download pre-built binaries from:
- https://sourceforge.net/projects/bowtie-bio/files/bowtie/

## Index Setup

ProbeDesign requires two types of bowtie indexes:
1. **Pseudogene indexes** - For `--pseudogene-mask` (e.g., `humanPseudo`)
2. **Genome indexes** - For `--genome-mask` (e.g., `GRCh38`)

### Pseudogene Indexes

Pseudogene FASTA files are included in `probedesign/pseudogeneDBs/`. Build indexes from them:

```bash
# Create indexes directory
mkdir -p bowtie_indexes
cd bowtie_indexes

# Build human pseudogene index
bowtie-build ../probedesign/pseudogeneDBs/human.fasta humanPseudo

# Build mouse pseudogene index (if needed)
bowtie-build ../probedesign/pseudogeneDBs/mouse.fasta mousePseudo

# Build other species as needed
bowtie-build ../probedesign/pseudogeneDBs/elegans.fasta celegansPseudo
bowtie-build ../probedesign/pseudogeneDBs/drosophila.fasta drosophilaPseudo
```

The index names must match what the code expects:
| Species | Index Name |
|---------|------------|
| human | `humanPseudo` |
| mouse | `mousePseudo` |
| elegans | `celegansPseudo` |
| drosophila | `drosophilaPseudo` |
| rat | `ratPseudo` |

### Human Genome (GRCh38)

1. **Create indexes directory**:
   ```bash
   mkdir -p bowtie_indexes
   cd bowtie_indexes
   ```

2. **Download the GRCh38 index** (~2.7GB compressed, ~3GB uncompressed):
   ```bash
   curl -O ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/GRCh38_no_alt.zip
   unzip GRCh38_no_alt.zip
   rm GRCh38_no_alt.zip  # Optional: remove zip to save space
   ```

3. **Verify the index files**:
   ```bash
   ls bowtie_indexes/
   # Should show (for human):
   # GCA_000001405.15_GRCh38_no_alt_analysis_set.*.ebwt  (genome)
   # humanPseudo.*.ebwt  (pseudogenes)
   ```

### Complete Setup for Human

After setup, your `bowtie_indexes/` directory should contain:
```
bowtie_indexes/
├── humanPseudo.1.ebwt          # Built from probedesign/pseudogeneDBs/human.fasta
├── humanPseudo.2.ebwt
├── humanPseudo.3.ebwt
├── humanPseudo.4.ebwt
├── humanPseudo.rev.1.ebwt
├── humanPseudo.rev.2.ebwt
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.1.ebwt   # Downloaded
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.2.ebwt
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.3.ebwt
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.4.ebwt
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.rev.1.ebwt
└── GCA_000001405.15_GRCh38_no_alt_analysis_set.rev.2.ebwt
```

### Other Available Indexes

Pre-built indexes are available at:
- https://bowtie-bio.sourceforge.net/index.shtml (click "Indexes")
- ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/

Common indexes:
- `GRCh38_no_alt.zip` - Human (GRCh38, recommended)
- `hg19.ebwt.zip` - Human (hg19, legacy)
- `mm10.ebwt.zip` - Mouse (mm10)
- `mm9.ebwt.zip` - Mouse (mm9, legacy)

## Environment Setup

For the ProbeDesign scripts to find bowtie, ensure it's in your PATH:

```bash
# If using conda/mamba, activate the environment:
conda activate base

# Verify bowtie is accessible:
which bowtie
```

The existing `DesignServer/bowtie_search.py` searches for bowtie in these locations:
1. `/usr/bin/bowtie`
2. `$BOWTIEHOME/bowtie` (set via environment variable)
3. `$HOME/bowtie/bowtie`
4. `$HOME/Dropbox (RajLab)/probeDesign/bowtie/bowtie`
5. `$HOME/Downloads/bowtie/bowtie`

You can set the `BOWTIEHOME` environment variable to point to your bowtie installation:
```bash
export BOWTIEHOME=/opt/homebrew/Caskroom/miniforge/base/bin
```

## Testing the Installation

```bash
# Test bowtie with a simple alignment
echo -e ">test\nACGTACGTACGTACGT" | bowtie -f -v 0 -k 1 GCA_000001405.15_GRCh38_no_alt_analysis_set -
```

## Disk Space Requirements

| Index | Compressed | Uncompressed |
|-------|------------|--------------|
| Human GRCh38 | ~2.7 GB | ~3.0 GB |
| Mouse mm10 | ~2.5 GB | ~2.7 GB |

## Troubleshooting

### "command not found: bowtie"
- Ensure conda environment is activated: `conda activate base`
- Check if bowtie is installed: `conda list | grep bowtie`

### "Could not find the bowtie executable"
- Set `BOWTIEHOME` environment variable
- Or install bowtie to one of the default search paths

### Index not found
- Ensure the index files are in the correct directory
- The index name is the prefix before `.1.ebwt` (e.g., `GCA_000001405.15_GRCh38_no_alt_analysis_set`)

## References

- Bowtie homepage: http://bowtie-bio.sourceforge.net/
- Bowtie manual: http://bowtie-bio.sourceforge.net/manual.shtml
- Bowtie GitHub: https://github.com/BenLangmead/bowtie
- Bioconda recipe: https://bioconda.github.io/recipes/bowtie/README.html
