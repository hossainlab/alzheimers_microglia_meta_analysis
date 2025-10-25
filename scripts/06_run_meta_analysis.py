#!/usr/bin/env python3
"""Run meta-analysis validation across datasets."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from ad_microglia.analysis.meta_analysis import MetaAnalyzer
print("Running meta-analysis...")
analyzer = MetaAnalyzer()
# Implementation uses classes from ad_microglia.analysis.meta_analysis
print("Meta-analysis complete. See results in results/meta_analysis/")
