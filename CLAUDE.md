# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the Raj Lab single molecule RNA FISH probe design software, primarily written in MATLAB with Python support scripts. The software designs specific oligonucleotide probes for targeting RNA sequences using various masking strategies to avoid off-target binding.

## Architecture

The codebase is organized into three main components:

### Core Probe Design (`probedesign/`)
- **Primary functions**: `findprobesHD.m` (server-based) and `findprobesLocal.m` (local bowtie)
- **Key algorithms**: Thermodynamic calculations, GC content analysis, repeat masking
- **Masking strategies**: Pseudogene, genome, repeat, GC run, and miRNA masking
- **Dependencies**: MATLAB, bowtie 0.12.7 for sequence alignment

### Design Server (`DesignServer/`)
- **Command-line interface**: `find_probes_cl.py` for batch probe design
- **Python modules**: `fasta.py`, `find_probes.py`, `probe_design.py`, `bowtie_search.py`
- **Bowtie integration**: Local alignment scripts and database management

### Pan-Probe Design (`panprobedesign/`)
- **Multi-species probes**: `findpanprobesHD.m` for designing probes across species
- **Alignment handling**: Functions for processing cross-species sequence alignments

## Common Development Tasks

### Running Probe Design
MATLAB commands for probe design:
```matlab
% Basic probe design (48 probes, default)
findprobesLocal('mysequence.fa')

% Custom number of probes
findprobesLocal('mysequence.fa', 32)

% Species-specific design with options
findprobesLocal('mysequence.fa', 32, 'species', 'mouse', 'targetTM', 72)
```

Python command-line interface:
```bash
python find_probes_cl.py input_sequence.fa human 48
```

### Bowtie Setup Requirements
- bowtie 0.12.7 executable must be in: `/usr/bin/bowtie`, `$HOME/bowtie`, `$HOME/Dropbox (RajLab)/probeDesign/bowtie`, or `$HOME/Downloads/bowtie`
- Indexed databases must be in `/indexes` subdirectory
- Supported species databases: human, mouse, elegans, drosophila (with genome, pseudogene, and mitochondrial variants)

### Python Environment
- Python 2.6+ required for DesignServer components
- Key modules: `fasta.py`, `find_probes.py`, `bowtie_search.py`
- No package manager files (requirements.txt, setup.py) - dependencies managed manually

## File Structure Notes

- **pseudogeneDBs/**: Contains species-specific pseudogene databases from pseudogene.org
- **src/**: Legacy web server dependencies (CherryPy, mod_python, soaplib)
- **maskprobes/**: Standalone masking utilities
- No build system - direct MATLAB/Python execution
- No automated tests - manual validation typical

## Important Functions

- `findprobesLocal.m`: Main local probe design function (replacement for server-based `findprobesHD.m`)
- `find_python_exe.m`: Cross-platform Python executable detection
- `hits_to_mask_local.m`: Local bowtie alignment for masking
- `pseudoTSVtoFASTA.m`: Convert pseudogene.org TSV files to FASTA format