# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2025-10-25

### Added
- Initial project structure with proper Python packaging
- Modular architecture with src/ layout
- Configuration management system
- Data preprocessing pipeline
- Cell type annotation using CellTypist
- Dataset integration (scVI and Scanorama)
- Activation state analysis
- Meta-analysis validation
- Comprehensive plot management system
- Test infrastructure
- Documentation framework

### Changed
- Restructured from flat scripts to proper package
- Centralized configuration in YAML files
- Improved import system (no more sys.path hacks)
- Better separation of concerns (src vs scripts)

### Infrastructure
- Added pyproject.toml for modern Python packaging
- Added setup.py for backward compatibility
- Comprehensive .gitignore
- Development requirements (pytest, black, mypy, etc.)
- Environment variable management (.env.example)
