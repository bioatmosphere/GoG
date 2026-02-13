#!/usr/bin/env python3
"""
Run GAPPY with DEMENTpy integration enabled.

This script demonstrates the complete integrated model with mechanistic
soil decomposition.
"""

import sys
from pathlib import Path

# Add model directory to path
sys.path.insert(0, str(Path(__file__).parent / 'model'))

from vegetation.gappy import GAPPYModel

def main():
    """Run GAPPY with DEMENTpy integration."""

    print("="*80)
    print("GAPPY-DEMENTpy Integrated Model")
    print("="*80)
    print()

    print("Initializing model...")
    model = GAPPYModel()

    # Enable DEMENTpy mechanistic soil decomposition
    print("Enabling DEMENTpy integration...")
    model.parameters.use_dement = True
    model.parameters.dement_spatial_mode = 'layered'

    # Reduce simulation length for testing
    print("Configuring simulation parameters...")
    model.parameters.numyears = 50  # Run 50 years for demonstration
    model.parameters.year_print_interval = 10
    model.parameters.spinup = False  # Skip spinup for faster testing

    print()
    print("Configuration:")
    print(f"  Simulation years: {model.parameters.numyears}")
    print(f"  Number of plots: {model.parameters.numplots}")
    print(f"  DEMENTpy enabled: {model.parameters.use_dement}")
    print(f"  Spatial mode: {model.parameters.dement_spatial_mode}")
    print(f"  Print interval: {model.parameters.year_print_interval} years")
    print()

    print("-"*80)
    print("Initializing input files and DEMENTpy grids...")
    print("-"*80)
    try:
        model.initialize_input_files()
        print("✓ Initialization complete")
        print()

        # Check that DEMENTpy was initialized
        if len(model.sites) > 0 and hasattr(model.sites[0].soil, 'dement_adapter'):
            adapter = model.sites[0].soil.dement_adapter
            if adapter:
                print("DEMENTpy Status:")
                print(f"  ✓ Adapter initialized")
                print(f"  ✓ Enabled: {adapter.enable_dement}")
                print(f"  ✓ Spatial mode: {adapter.spatial_mode}")
                print(f"  ✓ Number of grids: {adapter.n_grids}")
                if adapter.state_initialized:
                    print(f"  ✓ Grids initialized: {len(adapter.dement_grids)}")
                print()

    except Exception as e:
        print(f"✗ Initialization failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

    print("-"*80)
    print("Running simulation...")
    print("-"*80)
    try:
        # Run the model
        model.run()

        print()
        print("="*80)
        print("✓ Simulation completed successfully!")
        print("="*80)
        print()

        # Show coupling statistics if available
        if len(model.sites) > 0 and hasattr(model.sites[0].soil, 'dement_adapter'):
            adapter = model.sites[0].soil.dement_adapter
            if adapter and hasattr(adapter, 'coupling_stats'):
                stats = adapter.coupling_stats
                print("DEMENTpy Coupling Statistics:")
                print(f"  Total coupling calls: {stats['total_calls']}")
                print(f"  Total litter C: {stats['total_litter_c']:.2f} tc/ha")
                print(f"  Total litter N: {stats['total_litter_n']:.2f} tn/ha")
                print(f"  Total respiration: {stats['total_resp']:.2f} tc/ha")
                print(f"  Total N available: {stats['total_n_avail']:.2f} tn/ha")
                print()

        print("Output files written to: output_data/")
        print()

        return 0

    except Exception as e:
        print()
        print("="*80)
        print("✗ Simulation failed")
        print("="*80)
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
