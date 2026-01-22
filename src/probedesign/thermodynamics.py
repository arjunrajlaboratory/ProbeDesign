"""Thermodynamic calculations for RNA-DNA hybrids.

Uses Sugimoto et al. 1995 parameters for RNA-DNA hybridization and
SantaLucia PNAS 1998 Tm prediction equation with salt adjustment.

References:
- Sugimoto et al. (1995) "Thermodynamic Parameters To Predict Stability of
  RNA/DNA Hybrid Duplexes" Biochemistry 1995
- SantaLucia (1998) PNAS for Tm prediction
- OligoCalc website for salt adjustment
"""

import math
from typing import Tuple

# Sugimoto 1995 RNA-DNA parameters (Table 3)
# Enthalpy (kcal/mol) for dinucleotide pairs
# Index by first nucleotide (row) and second nucleotide (column)
# Format: dinucleotide 'XY' where X is 5' base, Y is 3' base
DELTA_H = {
    'aa': -7.8,  'ac': -5.9,  'ag': -9.1,  'at': -8.3,
    'ca': -9.0,  'cc': -9.3,  'cg': -16.3, 'ct': -7.0,
    'ga': -5.5,  'gc': -8.0,  'gg': -12.8, 'gt': -7.8,
    'ta': -7.8,  'tc': -8.6,  'tg': -10.4, 'tt': -11.5,
}

# Entropy (cal/(mol*K)) for dinucleotide pairs
DELTA_S = {
    'aa': -21.9, 'ac': -12.3, 'ag': -23.5, 'at': -23.9,
    'ca': -26.1, 'cc': -23.2, 'cg': -47.1, 'ct': -19.7,
    'ga': -13.5, 'gc': -17.1, 'gg': -31.9, 'gt': -21.6,
    'ta': -23.2, 'tc': -22.9, 'tg': -28.4, 'tt': -36.4,
}

# Initialization parameters
INIT_H = 1.9    # kcal/mol
INIT_S = -3.9   # cal/(mol*K)

# Standard conditions
DEFAULT_SALT = 0.33      # Sodium molarity (M) - 2X SSC buffer
DEFAULT_PRIMER_CONC = 50e-6  # Primer concentration (M) - 50 µM
R = 1.9872  # Gas constant in cal/(mol*K)


def thermo_rna_dna(seq: str) -> Tuple[float, float, float, float]:
    """Calculate thermodynamic values for RNA-DNA hybrid.

    Uses Sugimoto 1995 parameters for nearest-neighbor stacking.

    Args:
        seq: DNA sequence (5'->3'), will be converted to lowercase.
             Should contain only A, C, G, T characters.

    Returns:
        Tuple of (Tm, dG, dH, dS) where:
        - Tm: Melting temperature in °C
        - dG: Gibbs free energy in kcal/mol
        - dH: Enthalpy in kcal/mol
        - dS: Entropy in cal/(mol*K)
    """
    seq = seq.lower()

    # Sum enthalpy and entropy contributions from base stacking
    dH = 0.0  # kcal/mol
    dS = 0.0  # cal/(mol*K)

    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        dH += DELTA_H[dinuc]
        dS += DELTA_S[dinuc]

    # Add initialization terms
    dH += INIT_H
    dS += INIT_S

    # Calculate Gibbs free energy at 37°C (310.15 K)
    # dG = dH - T*dS, but dH is in kcal/mol and dS in cal/(mol*K)
    # So: dG = dH*1000 - T*dS (all in cal/mol), then convert back
    dG = (dH * 1000 - (37 + 273.15) * dS) / 1000  # kcal/mol

    # Calculate Tm using SantaLucia 1998 equation with salt adjustment
    # Tm = dH*1000 / (dS + R*ln(Ct/4)) - 273.15 + 16.6*log10(salt)
    Tm = (dH * 1000 / (dS + R * math.log(DEFAULT_PRIMER_CONC / 4))) - 273.15 + 16.6 * math.log10(DEFAULT_SALT)

    return Tm, dG, dH, dS


def gibbs_rna_dna(seq: str) -> float:
    """Calculate Gibbs free energy for RNA-DNA hybrid.

    Args:
        seq: DNA sequence (5'->3')

    Returns:
        Gibbs free energy in kcal/mol (negative values indicate favorable binding)
    """
    _, dG, _, _ = thermo_rna_dna(seq)
    return dG


def tm_rna_dna(seq: str) -> float:
    """Calculate melting temperature for RNA-DNA hybrid.

    Args:
        seq: DNA sequence (5'->3')

    Returns:
        Melting temperature in °C
    """
    Tm, _, _, _ = thermo_rna_dna(seq)
    return Tm


# SantaLucia 1998 DNA-DNA parameters (for optional DNA mode)
DNA_DELTA_H = {
    'aa': -7.9, 'ac': -8.4, 'ag': -7.8, 'at': -7.2,
    'ca': -8.5, 'cc': -8.0, 'cg': -10.6, 'ct': -7.8,
    'ga': -8.2, 'gc': -9.8, 'gg': -8.0, 'gt': -8.4,
    'ta': -7.2, 'tc': -8.2, 'tg': -8.5, 'tt': -7.9,
}

DNA_DELTA_S = {
    'aa': -22.2, 'ac': -22.4, 'ag': -21.0, 'at': -20.4,
    'ca': -22.7, 'cc': -19.9, 'cg': -27.2, 'ct': -21.0,
    'ga': -22.2, 'gc': -24.4, 'gg': -19.9, 'gt': -22.4,
    'ta': -21.3, 'tc': -22.2, 'tg': -22.7, 'tt': -22.2,
}


def tm_dna(seq: str) -> float:
    """Calculate melting temperature for DNA-DNA duplex.

    Uses SantaLucia 1998 parameters.

    Args:
        seq: DNA sequence (5'->3')

    Returns:
        Melting temperature in °C
    """
    seq = seq.lower()

    dH = 0.0
    dS = 0.0

    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        dH += DNA_DELTA_H[dinuc]
        dS += DNA_DELTA_S[dinuc]

    # Terminal corrections based on end nucleotides
    if seq[0] in ('c', 'g'):
        dH += 0.1
        dS += -2.8
    else:
        dH += 2.3
        dS += 4.1

    if seq[-1] in ('c', 'g'):
        dH += 0.1
        dS += -2.8
    else:
        dH += 2.3
        dS += 4.1

    Tm = (dH * 1000 / (dS + R * math.log(DEFAULT_PRIMER_CONC / 4))) + 16.6 * math.log10(DEFAULT_SALT) - 273.15

    return Tm
