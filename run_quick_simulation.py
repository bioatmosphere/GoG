#!/usr/bin/env python3
"""
Quick simulation run with DEMENTpy integration.
"""

import sys
import os

# Change to project directory
os.chdir('/Users/6lw/Desktop/2_models/GoG/GoG')

# Run using the same approach as test_integration.py
from model.vegetation.gappy import GAPPYModel

print("="*80)
print("GAPPY-DEMENTpy Integrated Model - Quick Simulation")
print("="*80)
print()

print("Initializing model with DEMENTpy...")
model = GAPPYModel()

# Enable DEMENTpy
model.parameters.use_dement = True
model.parameters.dement_spatial_mode = 'aggregated'

# Very short simulation for testing
model.parameters.numyears = 10
model.parameters.year_print_interval = 5
model.parameters.spinup = False

print(f"Configuration:")
print(f"  Years: {model.parameters.numyears}")
print(f"  Plots: {model.parameters.numplots}")
print(f"  DEMENTpy: {model.parameters.use_dement}")
print(f"  Mode: {model.parameters.dement_spatial_mode}")
print()

print("Initializing...")
model.initialize_input_files()

print()
print("✓ Model initialized")

# Check DEMENTpy status
if len(model.sites) > 0 and hasattr(model.sites[0].soil, 'dement_adapter'):
    adapter = model.sites[0].soil.dement_adapter
    if adapter:
        print(f"✓ DEMENTpy adapter present")
        print(f"  - Enabled: {adapter.enable_dement}")
        print(f"  - Grids: {adapter.n_grids}")
        print(f"  - Initialized: {adapter.state_initialized}")

print()
print("="*80)
print("Running 10-year simulation...")
print("="*80)

try:
    model.run()
    print()
    print("="*80)
    print("✓ SIMULATION COMPLETE!")
    print("="*80)

    # Show statistics
    if len(model.sites) > 0 and hasattr(model.sites[0].soil, 'dement_adapter'):
        adapter = model.sites[0].soil.dement_adapter
        if adapter:
            stats = adapter.coupling_stats
            print()
            print("Coupling Statistics:")
            print(f"  Calls: {stats['total_calls']}")
            print(f"  Litter C: {stats['total_litter_c']:.2f} tc/ha")
            print(f"  Litter N: {stats['total_litter_n']:.2f} tn/ha")
            print(f"  Respiration: {stats['total_resp']:.2f} tc/ha")
            print(f"  N available: {stats['total_n_avail']:.2f} tn/ha")

except Exception as e:
    print(f"\n✗ Simulation failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
