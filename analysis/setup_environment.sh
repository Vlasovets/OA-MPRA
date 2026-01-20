#!/bin/bash
# Setup script for MPRA analysis environment

echo "Creating conda environment for MPRA analysis..."

# Create the conda environment
conda env create -f environment.yml

echo ""
echo "Activating environment..."
conda activate mpra_analysis
