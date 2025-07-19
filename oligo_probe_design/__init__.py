"""
Oligo Probe Design Package

A Python package for designing specific oligonucleotide probes for single molecule RNA FISH.
This package provides a local, self-contained implementation of the probe design algorithm
originally developed in the Raj Lab.

Main Functions:
    design_probes: Main function for designing RNA FISH probes
"""

from .core.probe_designer import design_probes

__version__ = "1.0.0"
__author__ = "Raj Lab"

__all__ = ["design_probes"]