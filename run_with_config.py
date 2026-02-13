#!/usr/bin/env python3
"""
Run GAPPY with DEMENTpy using gappy_config.json.
"""

import sys
import os

# Change to project directory
os.chdir('/Users/6lw/Desktop/2_models/GoG/GoG')
sys.path.insert(0, os.path.join(os.getcwd(), 'model'))

from vegetation.gappy import GAPPYModel

print("="*80)
print("GAPPY-DEMENTpy Integrated Model")
print("="*80)
print()

print("Loading configuration from: input_data/gappy_config.json")
print()

# Create model (will load config automatically)
model = GAPPYModel()

# Initialize
print("Initializing model and input files...")
model.initialize_input_files()

print()
print("Configuration loaded:")
print(f"  Years: {model.parameters.numyears}")
print(f"  Plots: {model.parameters.numplots}")
print(f"  Spinup: {model.parameters.spinup}")
print(f"  Print interval: {model.parameters.year_print_interval}")
print(f"  DEMENTpy enabled: {model.parameters.use_dement}")
if model.parameters.use_dement:
    print(f"  DEMENTpy mode: {model.parameters.dement_spatial_mode}")

# Check DEMENTpy status
print()
if len(model.sites) > 0 and hasattr(model.sites[0].soil, 'dement_adapter'):
    adapter = model.sites[0].soil.dement_adapter
    if adapter:
        print("DEMENTpy Adapter Status:")
        print(f"  ✓ Initialized")
        print(f"  ✓ Enabled: {adapter.enable_dement}")
        print(f"  ✓ Spatial mode: {adapter.spatial_mode}")
        print(f"  ✓ Number of grids: {adapter.n_grids}")
        if adapter.state_initialized and len(adapter.dement_grids) > 0:
            print(f"  ✓ Grids created: {len(adapter.dement_grids)}")
            grid = adapter.dement_grids[0]
            if not isinstance(grid, dict):  # Real Grid, not placeholder
                print(f"  ✓ Grid type: Real DEMENTpy Grid")
                print(f"     - Substrates: {grid.Substrates.shape}")
                print(f"     - Microbes: {len(grid.Microbes)} cells")

print()
print("="*80)
print(f"Running {model.parameters.numyears}-year simulation...")
print("="*80)

try:
    model.run()

    print()
    print("="*80)
    print("✓ SIMULATION COMPLETE!")
    print("="*80)
    print()

    # Show coupling statistics
    if len(model.sites) > 0 and hasattr(model.sites[0].soil, 'dement_adapter'):
        adapter = model.sites[0].soil.dement_adapter
        if adapter and hasattr(adapter, 'coupling_stats'):
            stats = adapter.coupling_stats
            print("DEMENTpy Coupling Statistics:")
            print(f"  Total calls: {stats['total_calls']}")
            print(f"  Total litter C: {stats['total_litter_c']:.2f} tc/ha")
            print(f"  Total litter N: {stats['total_litter_n']:.4f} tn/ha")
            print(f"  Total respiration: {stats['total_resp']:.4f} tc/ha")
            print(f"  Total N available: {stats['total_n_avail']:.4f} tn/ha")
            print()

    print("Output files written to: output_data/")
    print()

except Exception as e:
    print()
    print("="*80)
    print("✗ SIMULATION FAILED")
    print("="*80)
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
