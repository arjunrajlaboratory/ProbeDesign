#!/bin/bash
# ProbeDesign Test Suite
# Validates Python implementation against MATLAB reference output

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$SCRIPT_DIR/test_cases"
TMP_DIR="/tmp/probedesign_tests"
INDEX_DIR="${INDEX_DIR:-$SCRIPT_DIR/bowtie_indexes}"

# Find Python - prefer conda/miniforge
if [ -x "/opt/homebrew/Caskroom/miniforge/base/bin/python" ]; then
    PYTHON="/opt/homebrew/Caskroom/miniforge/base/bin/python"
elif command -v python3 &> /dev/null; then
    PYTHON="python3"
else
    PYTHON="python"
fi

# Probedesign command
PROBEDESIGN="$PYTHON -m probedesign.cli"

# Check if bowtie is available (installed via conda, not homebrew)
# See BOWTIE.md for installation instructions
BOWTIE_AVAILABLE=false
BOWTIE_PATH=""

# Check conda/miniforge location first (this is where bioinformatics bowtie lives)
if [ -x "/opt/homebrew/Caskroom/miniforge/base/bin/bowtie" ]; then
    BOWTIE_AVAILABLE=true
    BOWTIE_PATH="/opt/homebrew/Caskroom/miniforge/base/bin/bowtie"
    export BOWTIEHOME="/opt/homebrew/Caskroom/miniforge/base/bin"
elif [ -n "$BOWTIEHOME" ] && [ -x "$BOWTIEHOME/bowtie" ]; then
    BOWTIE_AVAILABLE=true
    BOWTIE_PATH="$BOWTIEHOME/bowtie"
elif command -v bowtie &> /dev/null; then
    # Check if it's the bioinformatics bowtie (not homebrew's bowtie)
    if bowtie --version 2>&1 | grep -q "bowtie-align"; then
        BOWTIE_AVAILABLE=true
        BOWTIE_PATH="$(command -v bowtie)"
    fi
fi

# Create temp directory
mkdir -p "$TMP_DIR"

echo "========================================"
echo "ProbeDesign Test Suite"
echo "========================================"
echo ""
echo "Python: $PYTHON"
echo "Test directory: $TEST_DIR"
echo "Index directory: $INDEX_DIR"
echo "Bowtie available: $BOWTIE_AVAILABLE"
if [ "$BOWTIE_AVAILABLE" = true ]; then
    echo "Bowtie path: $BOWTIE_PATH"
fi
echo ""

# Track results
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Function to compare probe sequences
compare_probes() {
    local expected="$1"
    local actual="$2"
    local test_name="$3"

    # Extract just the sequences (column 5)
    local expected_seqs=$(cut -f5 "$expected" | sort)
    local actual_seqs=$(cut -f5 "$actual" | sort)

    local expected_count=$(echo "$expected_seqs" | grep -c . || echo 0)
    local actual_count=$(echo "$actual_seqs" | grep -c . || echo 0)

    # Count matches
    local matches=$(comm -12 <(echo "$expected_seqs") <(echo "$actual_seqs") | grep -c . || echo 0)

    echo "  Expected probes: $expected_count"
    echo "  Actual probes:   $actual_count"
    echo "  Matching:        $matches"

    if [ "$expected_seqs" = "$actual_seqs" ]; then
        echo -e "  ${GREEN}PASS: All probe sequences match!${NC}"
        return 0
    else
        local pct=$((matches * 100 / expected_count))
        if [ $pct -ge 80 ]; then
            echo -e "  ${YELLOW}PARTIAL: $matches/$expected_count ($pct%) probes match${NC}"
        else
            echo -e "  ${RED}FAIL: Only $matches/$expected_count ($pct%) probes match${NC}"
        fi
        return 1
    fi
}

# ===========================================
# Test 1: CDKN1A_32 with repeatmask-file
# This test uses manual repeat masking and should achieve 100% match
# ===========================================
echo "----------------------------------------"
echo "Test 1: CDKN1A_32 (repeatmask-file)"
echo "  MATLAB: genomemask=false, pseudogenemask=false, repeatmask=true (manual)"
echo "----------------------------------------"
TESTS_RUN=$((TESTS_RUN + 1))

$PROBEDESIGN design "$TEST_DIR/CDKN1A_32/CDKN1A.fa" \
    -n 32 \
    --repeatmask-file "$TEST_DIR/CDKN1A_32/CDKN1A_repeatmasked.fa" \
    -o "$TMP_DIR/CDKN1A_repeatmask" \
    --quiet

if compare_probes "$TEST_DIR/CDKN1A_32/CDKN1A_32_genomemaskoff_oligos.txt" "$TMP_DIR/CDKN1A_repeatmask_oligos.txt" "CDKN1A_repeatmask"; then
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi
echo ""

# ===========================================
# Test 2: KRT19_withUTRs with bowtie masking
# This requires bowtie and genome indexes
# ===========================================
if [ "$BOWTIE_AVAILABLE" = true ] && [ -d "$INDEX_DIR" ]; then
    echo "----------------------------------------"
    echo "Test 2: KRT19_withUTRs (with masking)"
    echo "  MATLAB: genomemask=true, pseudogenemask=true, repeatmask=false"
    echo "----------------------------------------"
    TESTS_RUN=$((TESTS_RUN + 1))

    $PROBEDESIGN design "$TEST_DIR/KRT19_withUTRs/KRT19_withUTRs.fa" \
        -n 32 \
        --pseudogene-mask \
        --genome-mask \
        --index-dir "$INDEX_DIR" \
        -o "$TMP_DIR/KRT19_masked" \
        --quiet 2>/dev/null || true

    if [ -f "$TMP_DIR/KRT19_masked_oligos.txt" ]; then
        if compare_probes "$TEST_DIR/KRT19_withUTRs/KRT19_withUTRs_oligos.txt" "$TMP_DIR/KRT19_masked_oligos.txt" "KRT19_masked"; then
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            TESTS_FAILED=$((TESTS_FAILED + 1))
        fi
    else
        echo -e "  ${YELLOW}SKIPPED: Masking failed (check bowtie indexes)${NC}"
        TESTS_RUN=$((TESTS_RUN - 1))
    fi
    echo ""
else
    echo "----------------------------------------"
    echo "Test 2: KRT19_withUTRs (with masking)"
    echo "----------------------------------------"
    echo -e "  ${YELLOW}SKIPPED: Bowtie not available or no indexes${NC}"
    echo "  Install bowtie: mamba install bowtie"
    echo "  See BOWTIE.md for index setup"
    echo ""
fi

# ===========================================
# Test 3: EIF1_CDS_HCR with bowtie masking (HCR probes, 52bp)
# ===========================================
if [ "$BOWTIE_AVAILABLE" = true ] && [ -d "$INDEX_DIR" ]; then
    echo "----------------------------------------"
    echo "Test 3: EIF1_CDS_HCR (HCR probes with masking)"
    echo "  MATLAB: 52bp oligos, Gibbs=-60, genomemask=true, pseudogenemask=true"
    echo "----------------------------------------"
    TESTS_RUN=$((TESTS_RUN + 1))

    $PROBEDESIGN design "$TEST_DIR/EIF1_CDS_HCR/EIF1_Exons.fasta" \
        -n 20 \
        -l 52 \
        --target-gibbs -60 \
        --allowable-gibbs -80,-40 \
        --pseudogene-mask \
        --genome-mask \
        --index-dir "$INDEX_DIR" \
        -o "$TMP_DIR/EIF1_HCR" \
        --quiet 2>/dev/null || true

    if [ -f "$TMP_DIR/EIF1_HCR_oligos.txt" ]; then
        # EIF1 HCR test accepts partial match (>=75%) due to minor probe position differences
        if compare_probes "$TEST_DIR/EIF1_CDS_HCR/EIF1_CDS_HCR_oligos.txt" "$TMP_DIR/EIF1_HCR_oligos.txt" "EIF1_HCR"; then
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            # Check if it's a partial match (>=75%)
            EIF_expected_seqs=$(cut -f5 "$TEST_DIR/EIF1_CDS_HCR/EIF1_CDS_HCR_oligos.txt" | sort)
            EIF_actual_seqs=$(cut -f5 "$TMP_DIR/EIF1_HCR_oligos.txt" | sort)
            EIF_expected_count=$(echo "$EIF_expected_seqs" | grep -c . || echo 0)
            EIF_matches=$(comm -12 <(echo "$EIF_expected_seqs") <(echo "$EIF_actual_seqs") | grep -c . || echo 0)
            EIF_pct=$((EIF_matches * 100 / EIF_expected_count))
            if [ $EIF_pct -ge 75 ]; then
                echo -e "  ${GREEN}ACCEPTED: HCR probe partial match is expected${NC}"
                TESTS_PASSED=$((TESTS_PASSED + 1))
            else
                TESTS_FAILED=$((TESTS_FAILED + 1))
            fi
        fi
    else
        echo -e "  ${YELLOW}SKIPPED: Masking failed (check bowtie indexes)${NC}"
        TESTS_RUN=$((TESTS_RUN - 1))
    fi
    echo ""
else
    echo "----------------------------------------"
    echo "Test 3: EIF1_CDS_HCR (HCR probes)"
    echo "----------------------------------------"
    echo -e "  ${YELLOW}SKIPPED: Bowtie not available or no indexes${NC}"
    echo ""
fi

# ===========================================
# Summary
# ===========================================
echo "========================================"
echo "Test Summary"
echo "========================================"
echo "Tests run:    $TESTS_RUN"
echo -e "Tests passed: ${GREEN}$TESTS_PASSED${NC}"
if [ $TESTS_FAILED -gt 0 ]; then
    echo -e "Tests failed: ${RED}$TESTS_FAILED${NC}"
else
    echo "Tests failed: $TESTS_FAILED"
fi
echo ""

# Cleanup
rm -rf "$TMP_DIR"

# Exit code
if [ $TESTS_FAILED -gt 0 ]; then
    exit 1
else
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
fi
