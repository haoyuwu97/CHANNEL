"""
CHANNEL Python interface.

This module mirrors the ergonomics of PILOTS-style pipelines:
- run a simulation from a config.ini
- optionally wire in PILOTS kernel outputs
- load results in a structured way
"""

from .api import run, run_from_pilots, load_results, discover_pilots_kernels

__all__ = ["run", "run_from_pilots", "load_results", "discover_pilots_kernels"]
