"""
Pytest configuration and fixtures for AD microglia tests.
"""

import pytest
from pathlib import Path
import sys

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


@pytest.fixture
def project_root():
    """Return project root directory."""
    return Path(__file__).parent.parent


@pytest.fixture
def test_data_dir():
    """Return test data directory."""
    return Path(__file__).parent / "test_data"


@pytest.fixture
def sample_settings():
    """Return sample settings for testing."""
    from ad_microglia.config.settings import Settings
    return Settings()
