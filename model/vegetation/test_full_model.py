"""
Full model test for the complete UVAFME Python implementation.
"""

import sys
import os
sys.path.append('..')

from vegetation import UVAFMEModel, params

def test_full_model():
    """Test the complete UVAFME model with all components."""
    
    print("=" * 60)
    print("UVAFME Python Implementation - Full Model Test")
    print("=" * 60)
    print()
    
    # Create model instance
    model = UVAFMEModel()
    
    # Set up small test parameters
    params.numyears = 20
    params.numplots = 1
    params.maxtrees = 100
    params.maxheight = 40
    params.year_print_interval = 5
    params.fixed_seed = True
    params.same_climate = True
    
    print("Test configuration:")
    print(f"  Years: {params.numyears}")
    print(f"  Plots per site: {params.numplots}")
    print(f"  Max trees per plot: {params.maxtrees}")
    print(f"  Print interval: {params.year_print_interval}")
    print()
    
    try:
        # Run the model
        model.run()
        
        print()
        print("=" * 60)
        print("Full model test completed successfully! ✓")
        print("=" * 60)
        
        # Show output files created
        output_dir = "output_data"
        if os.path.exists(output_dir):
            print(f"\nOutput files created in {output_dir}:")
            for file in os.listdir(output_dir):
                filepath = os.path.join(output_dir, file)
                if os.path.isfile(filepath):
                    size = os.path.getsize(filepath)
                    print(f"  {file} ({size} bytes)")
        
    except Exception as e:
        print(f"Full model test failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True

if __name__ == "__main__":
    success = test_full_model()
    sys.exit(0 if success else 1)