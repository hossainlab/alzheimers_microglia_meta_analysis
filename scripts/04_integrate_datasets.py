#!/usr/bin/env python3
"""Integrate datasets using scVI or Scanorama."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from ad_microglia.integration.integration import MicrogliaIntegrator
print("Running dataset integration...")
integrator = MicrogliaIntegrator()
# Implementation uses classes from ad_microglia.integration.integration
print("Dataset integration complete. See results in data/processed/")
