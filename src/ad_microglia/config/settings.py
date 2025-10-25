"""
Centralized configuration settings for AD microglia meta-analysis.
"""

from pathlib import Path
from typing import Optional
import os


class Settings:
    """Global settings for the analysis pipeline."""

    def __init__(self, base_dir: Optional[Path] = None):
        """
        Initialize settings.

        Parameters
        ----------
        base_dir : Path, optional
            Base directory for the project. Defaults to project root.
        """
        if base_dir is None:
            # Auto-detect project root (where pyproject.toml or setup.py exists)
            current = Path(__file__).resolve()
            for parent in current.parents:
                if (parent / "pyproject.toml").exists() or (parent / "setup.py").exists():
                    base_dir = parent
                    break
            else:
                # Fallback to 3 levels up from this file
                base_dir = Path(__file__).resolve().parents[3]

        self.BASE_DIR = Path(base_dir)

        # Data directories
        self.DATA_DIR = self.BASE_DIR / "data"
        self.RAW_DATA_DIR = self.DATA_DIR / "raw"
        self.PROCESSED_DATA_DIR = self.DATA_DIR / "processed"
        self.METADATA_DIR = self.DATA_DIR / "metadata"

        # Results directories
        self.RESULTS_DIR = self.BASE_DIR / "results"
        self.FIGURES_DIR = self.RESULTS_DIR / "figures"
        self.TABLES_DIR = self.RESULTS_DIR / "tables"
        self.REPORTS_DIR = self.RESULTS_DIR / "reports"
        self.META_ANALYSIS_DIR = self.RESULTS_DIR / "meta_analysis"

        # Plot directories (legacy support)
        self.PLOTS_DIR = self.BASE_DIR / "plots"

        # Configuration directory
        self.CONFIG_DIR = self.BASE_DIR / "config"

        # Notebooks directory
        self.NOTEBOOKS_DIR = self.BASE_DIR / "notebooks"

        # Ensure critical directories exist
        self._create_directories()

        # Analysis parameters
        self.QC_MIN_GENES = int(os.getenv("QC_MIN_GENES", "200"))
        self.QC_MIN_CELLS = int(os.getenv("QC_MIN_CELLS", "3"))
        self.QC_MAX_MT_PERCENT = float(os.getenv("QC_MAX_MT_PERCENT", "20"))
        self.QC_MAX_GENES = int(os.getenv("QC_MAX_GENES", "5000"))
        self.QC_MIN_COUNTS = int(os.getenv("QC_MIN_COUNTS", "1000"))

        # Normalization parameters
        self.NORM_TARGET_SUM = float(os.getenv("NORM_TARGET_SUM", "10000"))

        # Integration parameters
        self.INTEGRATION_METHOD = os.getenv("INTEGRATION_METHOD", "scvi")

        # Plotting parameters
        self.PLOT_DPI = int(os.getenv("PLOT_DPI", "300"))
        self.PLOT_FORMAT = os.getenv("PLOT_FORMAT", "pdf,png").split(",")

        # Logging
        self.LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO")

        # Datasets
        self.DATASETS = {
            "GSE98969": "Keren-Shaul et al. (2017) - Mouse AD microglia",
            "GSE103334": "Mathys et al. (2019) - Human AD hippocampus",
            "GSE135437": "Sankowski et al. (2019) - Human AD cortex",
            "GSE157827": "Leng et al. (2021) - Human AD prefrontal cortex",
            "GSE129788": "Aging mouse brain microglia",
        }

    def _create_directories(self):
        """Create necessary directories if they don't exist."""
        directories = [
            self.DATA_DIR,
            self.RAW_DATA_DIR,
            self.PROCESSED_DATA_DIR,
            self.METADATA_DIR,
            self.RESULTS_DIR,
            self.FIGURES_DIR,
            self.TABLES_DIR,
            self.REPORTS_DIR,
            self.META_ANALYSIS_DIR,
        ]

        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)

    def get_raw_data_path(self, dataset_name: str) -> Path:
        """Get path to raw data for a dataset."""
        return self.RAW_DATA_DIR / dataset_name

    def get_processed_data_path(self, dataset_name: str, suffix: str = "processed") -> Path:
        """Get path to processed data file."""
        return self.PROCESSED_DATA_DIR / f"{dataset_name}_{suffix}.h5ad"

    def get_integrated_data_path(self, method: str = "scvi") -> Path:
        """Get path to integrated data file."""
        return self.PROCESSED_DATA_DIR / f"microglia_integrated_{method}.h5ad"


# Global settings instance
_settings: Optional[Settings] = None


def get_settings(base_dir: Optional[Path] = None) -> Settings:
    """
    Get global settings instance (singleton pattern).

    Parameters
    ----------
    base_dir : Path, optional
        Base directory for the project. Only used on first call.

    Returns
    -------
    Settings
        Global settings instance
    """
    global _settings
    if _settings is None:
        _settings = Settings(base_dir)
    return _settings
