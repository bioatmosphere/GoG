"""
Basic test script for the Python UVAFME vegetation model.
"""

import sys
sys.path.append('..')

from vegetation import SpeciesData, SiteData, SoilData, UVAFMEModel
from vegetation.constants import CODENAME, VERSION_ID

def test_species_data():
    """Test SpeciesData class."""
    print("Testing SpeciesData...")
    
    species = SpeciesData()
    
    # Test initialization
    species.initialize_species(
        species_id=1,
        genus_name="Quercus",
        taxonomic_name="Quercus alba",
        unique_id="QUAL",
        common_name="White Oak",
        genus_id=1,
        shade_tol=3,
        lownutr_tol=2,
        stress_tol=2,
        age_tol=5,
        drought_tol=2,
        flood_tol=2,
        fire_tol=2,
        max_age=300.0,
        max_diam=150.0,
        max_ht=35.0,
        wood_bulk_dens=0.68,
        rootdepth=2.0,
        leafdiam_a=0.1,
        leafarea_c=50.0,
        deg_day_min=800.0,
        deg_day_opt=1200.0,
        deg_day_max=1800.0,
        seedling_lg=0.5,
        invader=0.0,
        seed_num=100.0,
        sprout_num=50.0,
        seed_surv=0.1,
        arfa_0=0.8,
        g=0.3,
        conifer=False
    )
    
    # Test response functions
    light_response = species.light_rsp(0.5)
    species.temp_rsp(1200.0)
    species.drought_rsp(0.1, 0.05)
    species.flood_rsp(0.1)
    species.fire_rsp(0)
    
    print(f"  Species: {species.common_name}")
    print(f"  Light response (0.5): {light_response:.3f}")
    print(f"  Temperature response factor: {species.fc_degday:.3f}")
    print(f"  Drought response factor: {species.fc_drought:.3f}")
    print("  ✓ SpeciesData test passed")
    print()


def test_site_data():
    """Test SiteData class."""
    print("Testing SiteData...")
    
    site = SiteData()
    
    # Test basic initialization
    site.initialize_site(
        siteid=1,
        sitename="Test Site",
        siteregion="Virginia",
        lat=37.0,
        long=-78.0,
        wmo=123456,
        elevation=200.0,
        slope=5.0,
        Afc=20.0,
        A_perm_wp=8.0,
        lai=3.0,
        base_h=70.0,
        lai_w0=0.5,
        A0_w0=2.0,
        A_w0=15.0,
        sbase_w0=40.0,
        fire_prob=0.01,
        wind_prob=0.05,
        A0_c0=10.0,
        A0_n0=1.0,
        A_c0=50.0,
        A_n0=5.0,
        sbase_c0=100.0,
        sbase_n0=10.0,
        sigma=0.1,
        temp_lapse=[0.6] * 12,
        prcp_lapse=[0.1] * 12
    )
    
    # Test climate attachment
    tmin = [-5, -2, 3, 8, 13, 18, 20, 19, 15, 9, 4, -1]
    tmax = [5, 8, 15, 22, 27, 32, 35, 33, 28, 22, 14, 7]
    prcp = [80, 70, 90, 100, 120, 110, 90, 85, 75, 65, 70, 85]
    
    site.attach_climate(tmin, tmax, prcp)
    
    print(f"  Site: {site.site_name}")
    print(f"  Location: {site.latitude:.1f}°N, {site.longitude:.1f}°W")
    print(f"  Elevation: {site.elevation:.0f}m")
    print(f"  Annual temp range: {min(site.tmin):.1f}°C to {max(site.tmax):.1f}°C")
    print("  ✓ SiteData test passed")
    print()


def test_soil_data():
    """Test SoilData class."""
    print("Testing SoilData...")
    
    soil = SoilData()
    
    # Initialize with some test values
    soil.A0_c0 = 10.0
    soil.A0_n0 = 1.0
    soil.A_c0 = 50.0
    soil.A_n0 = 5.0
    soil.BL_c0 = 100.0
    soil.BL_n0 = 10.0
    soil.A_field_cap = 20.0
    soil.A_perm_wp = 8.0
    soil.A0_w0 = 2.0
    soil.A_w0 = 15.0
    soil.BL_w0 = 40.0
    
    # Test decomposition
    avail_N, C_resp = soil.soil_decomp(
        litter_c1=1.0,
        litter_c2=0.5,
        litter_n1=0.1,
        litter_n2=0.05,
        tempC=15.0,
        precip=2.0,
        aow0_scaled_by_max=0.5,
        saw0_scaled_by_fc=0.7,
        sbw0_scaled_by_max=0.8
    )
    
    print(f"  Available N: {avail_N:.4f}")
    print(f"  C respiration: {C_resp:.4f}")
    print("  ✓ SoilData test passed")
    print()


def test_model_creation():
    """Test basic model creation."""
    print("Testing UVAFMEModel...")
    
    model = UVAFMEModel()
    
    print(f"  Model type: {type(model).__name__}")
    print(f"  Code name: {CODENAME}")
    print(f"  Version: {VERSION_ID}")
    print("  ✓ UVAFMEModel creation test passed")
    print()


def main():
    """Run all tests."""
    print("=" * 60)
    print(f"UVAFME Python Translation - Basic Tests")
    print("=" * 60)
    print()
    
    try:
        test_species_data()
        test_site_data()
        test_soil_data()
        test_model_creation()
        
        print("=" * 60)
        print("All basic tests passed! ✓")
        print("=" * 60)
        
    except Exception as e:
        print(f"Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()