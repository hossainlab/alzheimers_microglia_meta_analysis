#!/usr/bin/env python3
"""Analyze microglial activation states."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from ad_microglia.analysis.activation_states import ActivationStateAnalyzer
print("Running activation state analysis...")
analyzer = ActivationStateAnalyzer()
# Implementation uses classes from ad_microglia.analysis.activation_states
print("Activation state analysis complete. See results in results/")
