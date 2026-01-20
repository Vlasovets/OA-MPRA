#!/bin/bash
# Setup script for MPRA analysis environment

echo "Creating conda environment for MPRA analysis..."

# Create the conda environment
conda env create -f environment.yml

echo ""
echo "Environment created successfully!"
echo "To activate the environment, run:"
echo "  conda activate mpra_analysis"
echo ""
echo "Then you can run the analysis scripts:"
echo "  cd /home/itg/oleg.vlasovets/projects/MPRA_data/mpra_test/analysis"
echo "  python 01_load_data.py"
echo "  python 02_mpra_analysis.py"
echo "  python 03_visualization.py"
