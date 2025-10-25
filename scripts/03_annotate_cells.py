#!/usr/bin/env python3
"""Annotate cell types to identify microglial cells."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from ad_microglia.analysis.cell_annotation import CellAnnotator
print("Running cell annotation...")
annotator = CellAnnotator()
# Implementation uses classes from ad_microglia.analysis.cell_annotation
print("Cell annotation complete. See results in data/processed/")
