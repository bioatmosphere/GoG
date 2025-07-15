#!/usr/bin/env python3
"""
Utility to create default UVAFME configuration file.
"""

import sys
import os
sys.path.append('.')

from vegetation import params

def create_default_config():
    """Create a default UVAFME configuration file."""
    
    # Create config in current directory
    config_file = "uvafme_config.json"
    params.create_default_config(config_file)
    
    # Also create in input_data directory
    input_dir = "input_data"
    os.makedirs(input_dir, exist_ok=True)
    input_config_file = os.path.join(input_dir, config_file)
    params.create_default_config(input_config_file)
    
    print(f"Created configuration files:")
    print(f"  - {config_file}")
    print(f"  - {input_config_file}")
    
    print("\nConfiguration parameters:")
    print(f"  Years to simulate: {params.numyears}")
    print(f"  Plots per site: {params.numplots}")
    print(f"  Max trees per plot: {params.maxtrees}")
    print(f"  Print interval: {params.year_print_interval} years")
    print(f"  Fixed seed: {params.fixed_seed}")
    print(f"  Plot size: {params.plotsize} ha")
    print(f"  Root depth: {params.rootdepth} m")
    
    print("\nYou can now run the model with: uv run python main.py")

if __name__ == "__main__":
    create_default_config()