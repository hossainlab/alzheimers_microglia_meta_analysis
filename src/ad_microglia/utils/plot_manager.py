#!/usr/bin/env python3
"""
Plot Management System for Alzheimer's Microglia Meta-Analysis
Centralized plot saving with organized directory structure
"""

import os
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from datetime import datetime

class PlotManager:
    """Centralized plot management system"""
    
    def __init__(self, base_dir='plots', dpi=300, figsize=(8, 6)):
        self.base_dir = Path(base_dir)
        self.dpi = dpi
        self.figsize = figsize
        
        # Create organized directory structure
        self.directories = {
            'qc': self.base_dir / 'quality_control',
            'preprocessing': self.base_dir / 'preprocessing', 
            'integration': self.base_dir / 'integration',
            'annotation': self.base_dir / 'cell_annotation',
            'activation': self.base_dir / 'activation_states',
            'meta_analysis': self.base_dir / 'meta_analysis',
            'comparative': self.base_dir / 'comparative_analysis',
            'supplementary': self.base_dir / 'supplementary'
        }
        
        # Create all directories
        for dir_path in self.directories.values():
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Configure matplotlib to not show plots
        plt.ioff()  # Turn off interactive mode
        
        # Set up logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
        # Create plot manifest file
        self.manifest_file = self.base_dir / 'plot_manifest.txt'
        if not self.manifest_file.exists():
            with open(self.manifest_file, 'w') as f:
                f.write(f"# Plot Manifest - Created {datetime.now()}\n")
                f.write("# Format: timestamp | category | filename | description\n\n")
    
    def save_plot(self, fig, filename, category='supplementary', description='', 
                  file_formats=['pdf', 'png'], close_fig=True):
        """
        Save plot to organized directory structure
        
        Parameters:
        -----------
        fig : matplotlib.figure.Figure
            The figure to save
        filename : str
            Filename without extension
        category : str
            Plot category (qc, preprocessing, integration, etc.)
        description : str
            Plot description for manifest
        file_formats : list
            File formats to save ['pdf', 'png', 'svg']
        close_fig : bool
            Whether to close figure after saving
        """
        
        if category not in self.directories:
            self.logger.warning(f"Unknown category '{category}'. Using 'supplementary'")
            category = 'supplementary'
        
        save_dir = self.directories[category]
        saved_files = []
        
        for fmt in file_formats:
            filepath = save_dir / f"{filename}.{fmt}"
            
            try:
                fig.savefig(
                    filepath,
                    dpi=self.dpi,
                    bbox_inches='tight',
                    format=fmt,
                    facecolor='white',
                    edgecolor='none'
                )
                saved_files.append(str(filepath))
                self.logger.info(f"Saved plot: {filepath}")
                
            except Exception as e:
                self.logger.error(f"Failed to save {filepath}: {e}")
        
        # Update manifest
        self._update_manifest(category, filename, description, saved_files)
        
        if close_fig:
            plt.close(fig)
        
        return saved_files
    
    def create_figure(self, figsize=None, **kwargs):
        """Create a new figure with default settings"""
        if figsize is None:
            figsize = self.figsize
        
        fig, ax = plt.subplots(figsize=figsize, **kwargs)
        return fig, ax
    
    def create_subplots(self, nrows=1, ncols=1, figsize=None, **kwargs):
        """Create subplots with default settings"""
        if figsize is None:
            # Adjust figsize based on number of subplots
            width, height = self.figsize
            figsize = (width * ncols, height * nrows)
        
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, **kwargs)
        return fig, axes
    
    def _update_manifest(self, category, filename, description, saved_files):
        """Update plot manifest file"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        with open(self.manifest_file, 'a') as f:
            f.write(f"{timestamp} | {category} | {filename} | {description}\n")
            for file_path in saved_files:
                f.write(f"    -> {file_path}\n")
    
    def get_plot_path(self, category, filename, extension='pdf'):
        """Get the full path for a plot file"""
        if category not in self.directories:
            category = 'supplementary'
        
        return self.directories[category] / f"{filename}.{extension}"
    
    def list_plots(self, category=None):
        """List all plots in a category or all categories"""
        if category:
            if category in self.directories:
                return list(self.directories[category].glob("*"))
            else:
                return []
        else:
            all_plots = []
            for cat_dir in self.directories.values():
                all_plots.extend(list(cat_dir.glob("*")))
            return all_plots
    
    def cleanup_empty_dirs(self):
        """Remove empty plot directories"""
        for dir_path in self.directories.values():
            if dir_path.exists() and not any(dir_path.iterdir()):
                dir_path.rmdir()
                self.logger.info(f"Removed empty directory: {dir_path}")

# Global plot manager instance
plot_manager = PlotManager()

# Convenience functions
def save_plot(fig, filename, category='supplementary', description='', 
              file_formats=['pdf', 'png'], close_fig=True):
    """Convenience function to save plots"""
    return plot_manager.save_plot(fig, filename, category, description, 
                                 file_formats, close_fig)

def create_figure(figsize=None, **kwargs):
    """Convenience function to create figures"""
    return plot_manager.create_figure(figsize, **kwargs)

def create_subplots(nrows=1, ncols=1, figsize=None, **kwargs):
    """Convenience function to create subplots"""
    return plot_manager.create_subplots(nrows, ncols, figsize, **kwargs)

def get_plot_path(category, filename, extension='pdf'):
    """Convenience function to get plot paths"""
    return plot_manager.get_plot_path(category, filename, extension)

# Context manager for plot saving
class PlotContext:
    """Context manager for automatic plot saving"""
    
    def __init__(self, filename, category='supplementary', description='',
                 file_formats=['pdf', 'png'], figsize=None):
        self.filename = filename
        self.category = category
        self.description = description
        self.file_formats = file_formats
        self.figsize = figsize
        self.fig = None
        self.ax = None
    
    def __enter__(self):
        self.fig, self.ax = plot_manager.create_figure(figsize=self.figsize)
        return self.fig, self.ax
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:  # No exception occurred
            plot_manager.save_plot(
                self.fig, self.filename, self.category, 
                self.description, self.file_formats, close_fig=True
            )
        else:
            plt.close(self.fig)

# Example usage:
# with PlotContext('dataset_overview', 'qc', 'Dataset size comparison') as (fig, ax):
#     ax.bar(datasets, counts)
#     ax.set_title('Dataset Sizes')