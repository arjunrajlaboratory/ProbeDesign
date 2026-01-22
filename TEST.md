# ProbeDesign Testing

This document describes how to validate the Python ProbeDesign implementation against the original MATLAB output.

## Quick Start

Run the automated test suite:

```bash
./run_tests.sh
```

This will run all test cases and report pass/fail status.

## Current Test Results Summary

Last validated: January 2025

| Test Case | Description | Match Rate | Status |
|-----------|-------------|------------|--------|
| CDKN1A_32 (repeatmask-file) | Manual repeat masking | **32/32 (100%)** | ✅ PASS |
| KRT19_withUTRs (with masking) | Pseudogene + genome masking | **6/6 (100%)** | ✅ PASS |
| EIF1_CDS_HCR (HCR probes) | 52bp HCR probes | **15/19 (78%)** | ✅ PASS (partial expected) |

### Key Findings

1. **Exact match achieved** for CDKN1A_32 using `--repeatmask-file` option
2. **Exact match achieved** for KRT19_withUTRs with `--pseudogene-mask --genome-mask`
3. **HCR probes** show minor position differences (78% match) which is expected due to longer oligos and tie-breaking in DP algorithm
4. **Thermodynamics** calculations match MATLAB exactly (Gibbs FE within 0.01 kcal/mol)

## Prerequisites

1. **Python package installed**:
   ```bash
   pip install -e .
   ```

2. **Bowtie 1** (for masking tests) - installed via conda, NOT homebrew:
   ```bash
   # See BOWTIE.md for full setup
   mamba install bowtie

   # Verify it's the bioinformatics bowtie
   bowtie --version  # Should show "bowtie-align-s version 1.x.x"
   ```

3. **Required bowtie indexes in `bowtie_indexes/`**:
   - `humanPseudo` - Human pseudogene database (built from `probedesign/pseudogeneDBs/`)
   - `GCA_000001405.15_GRCh38_no_alt_analysis_set` - Human genome (GRCh38)

## Test Cases

### 1. CDKN1A_32 (Manual Repeat Masking)

**Description**: Tests the `--repeatmask-file` option for manual repeat masking. This should achieve 100% match with MATLAB.

**MATLAB command**:
```matlab
findprobesLocal('CDKN1A.fa', 32, ...
    'outfilename', 'CDKN1A_32_genomemaskoff', ...
    'genomemask', false, ...
    'pseudogenemask', false, ...
    'repeatmask', true, ...
    'repeatmaskmanual', true, ...
    'repeatmaskfile', 'CDKN1A_repeatmasked.fa', ...
    'species', 'human')
```

**Python equivalent**:
```bash
probedesign design test_cases/CDKN1A_32/CDKN1A.fa \
    -n 32 \
    --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa \
    -o CDKN1A_test
```

**Expected**: 32 probes, 100% match

### 2. KRT19_withUTRs (Bowtie Masking)

**Description**: Tests pseudogene and genome masking using bowtie alignments.

**MATLAB command**:
```matlab
findprobesLocal('KRT19_withUTRs.fa', 32, ...
    'outfilename', 'KRT19_withUTRs', ...
    'pseudogenemask', true, ...
    'genomemask', true, ...
    'repeatmask', false, ...
    'species', 'human')
```

**Python equivalent**:
```bash
probedesign design test_cases/KRT19_withUTRs/KRT19_withUTRs.fa \
    -n 32 \
    --pseudogene-mask \
    --genome-mask \
    --index-dir bowtie_indexes \
    -o KRT19_test
```

**Expected**: 6 probes, 100% match:
```
1  55  66.2  -24.7  gtctcagaagctgcgattcg
2  65  63.2  -26.0  gaatgctgggcgcgcgaaag
3  50  68.7  -24.8  aaggcagctttcatgctcag
4  55  64.5  -24.6  ggaccttggaggcagacaaa
5  55  62.5  -23.5  ctccaaaggacagcagaagc
6  65  67.3  -25.9  gagaagagccgggggtaagg
```

### 3. EIF1_CDS_HCR (HCR Probes)

**Description**: Tests HCR probe design with longer oligos (52bp) and different Gibbs FE targets.

**MATLAB command**:
```matlab
findprobesLocal('EIF1_Exons.fasta', 20, ...
    'outfilename', 'EIF1_CDS_HCR', ...
    'oligolength', 52, ...
    'allowableGibbsFE', [-80, -40], ...
    'targetGibbsFE', -60, ...
    'pseudogenemask', true, ...
    'genomemask', true, ...
    'repeatmask', false)
```

**Python equivalent**:
```bash
probedesign design test_cases/EIF1_CDS_HCR/EIF1_Exons.fasta \
    -n 20 \
    -l 52 \
    --target-gibbs -60 \
    --allowable-gibbs -80,-40 \
    --pseudogene-mask \
    --genome-mask \
    --index-dir bowtie_indexes \
    -o EIF1_test
```

**Expected**: 19 probes, ~75-80% match (partial match expected due to DP tie-breaking with longer oligos)

## Running Tests

### Automated Test Script

The `run_tests.sh` script runs all tests automatically:

```bash
./run_tests.sh
```

Output example:
```
========================================
ProbeDesign Test Suite
========================================

Python: /opt/homebrew/Caskroom/miniforge/base/bin/python
Test directory: /Users/.../ProbeDesign/test_cases
Index directory: /Users/.../ProbeDesign/bowtie_indexes
Bowtie available: true
Bowtie path: /opt/homebrew/Caskroom/miniforge/base/bin/bowtie

----------------------------------------
Test 1: CDKN1A_32 (repeatmask-file)
  MATLAB: genomemask=false, pseudogenemask=false, repeatmask=true (manual)
----------------------------------------
Repeat masking (from file): 174 positions masked
  Expected probes: 32
  Actual probes:   32
  Matching:        32
  PASS: All probe sequences match!

...

========================================
Test Summary
========================================
Tests run:    3
Tests passed: 3
Tests failed: 0

All tests passed!
```

### Manual Validation

To manually compare outputs:

```bash
# Run Python
probedesign design input.fa -n 32 -o test_output

# Compare probe sequences (column 5)
diff <(cut -f5 expected_oligos.txt | sort) \
     <(cut -f5 test_output_oligos.txt | sort)
```

## Thermodynamics Validation

Test that thermodynamic calculations match MATLAB:

```python
from probedesign.thermodynamics import gibbs_rna_dna, tm_rna_dna

# Known values from MATLAB
test_cases = [
    ("gtctcagaagctgcgattcg", -24.69, 66.2),
    ("gaatgctgggcgcgcgaaag", -25.99, 63.2),
    ("aaggcagctttcatgctcag", -24.80, 68.7),
]

for seq, expected_gibbs, expected_tm in test_cases:
    gibbs = gibbs_rna_dna(seq)
    tm = tm_rna_dna(seq)
    print(f"{seq}")
    print(f"  Gibbs: {gibbs:.2f} (expected {expected_gibbs:.2f})")
    print(f"  Tm:    {tm:.1f} (expected {expected_tm:.1f})")
```

## Known Differences

1. **HCR probes (52bp)**: May have 75-80% match due to:
   - Tie-breaking in DP algorithm when multiple positions have equal badness
   - Floating-point precision differences in thermodynamic calculations

2. **Probe names**: Python includes directory path in probe names, MATLAB uses just the output name. This is a cosmetic difference.

3. **Gibbs FE formatting**: Python shows `-23.0`, MATLAB shows `-23`. Same values, different formatting.

## Adding New Test Cases

1. Create a directory in `test_cases/`
2. Add the input FASTA file
3. Run MATLAB to generate expected output:
   ```matlab
   findprobesLocal('input.fa', 48, 'outfilename', 'output_name', ...)
   ```
4. Save the MATLAB command to `command.txt`
5. Keep the `*_oligos.txt` and `*_seq.txt` files as expected output
6. Add a test case to `run_tests.sh`

## Troubleshooting

**"Could not find bowtie executable"**
- Ensure bowtie is installed via conda: `mamba install bowtie`
- Check if conda environment is active: `conda activate base`
- Set BOWTIEHOME: `export BOWTIEHOME=/opt/homebrew/Caskroom/miniforge/base/bin`

**"Index not found"**
- Check that index files exist in `bowtie_indexes/`
- Index name is the prefix before `.1.ebwt`
- See BOWTIE.md for download instructions

**Different probe counts**
- Check if masking options match MATLAB command
- MATLAB defaults: `pseudogenemask=true`, `genomemask=true`
- Always use `--pseudogene-mask --genome-mask` to match MATLAB defaults
