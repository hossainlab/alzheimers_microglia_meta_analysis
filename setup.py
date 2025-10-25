"""
Setup script for ad-microglia package.
For compatibility with older build systems.
Modern builds should use pyproject.toml.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

setup(
    name="ad-microglia",
    version="0.1.0",
    author="AD Microglia Research Team",
    description="Meta-analysis toolkit for Alzheimer's disease microglia single-cell RNA-seq studies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/alzheimers_microglia_meta_analysis",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "pandas>=1.5.0",
        "numpy>=1.21.0",
        "scipy>=1.9.0",
        "matplotlib>=3.6.0",
        "seaborn>=0.12.0",
        "scanpy>=1.9.0",
        "anndata>=0.8.0",
        "scvi-tools>=1.0.0",
        "scanorama>=1.7.0",
        "celltypist>=1.6.0",
        "decoupler>=1.4.0",
        "statsmodels>=0.14.0",
        "scikit-learn>=1.1.0",
        "GEOparse>=2.0.3",
        "requests>=2.28.0",
        "tqdm>=4.64.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
            "black>=23.0",
            "flake8>=6.0",
            "mypy>=1.0",
            "ipython>=8.0",
            "jupyter>=1.0",
        ],
        "docs": [
            "sphinx>=5.0",
            "sphinx-rtd-theme>=1.0",
        ],
    },
    include_package_data=True,
)
