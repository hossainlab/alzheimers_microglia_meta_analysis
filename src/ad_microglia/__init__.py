"""
Alzheimer's Disease Microglia Meta-Analysis Package

A comprehensive toolkit for analyzing microglial activation states across
Alzheimer's disease single-cell RNA-sequencing studies.
"""

__version__ = "0.1.0"
__author__ = "AD Microglia Research Team"

# Import key components for easy access
from ad_microglia.utils.plot_manager import PlotManager, PlotContext

__all__ = [
    "PlotManager",
    "PlotContext",
]
