"""
UVAFME (UVA Forest Model Enhanced) main module.
Translated from UVAFME.f90
"""

import sys
import time
import numpy as np
from .constants import *
from .parameters import params
from .species import SpeciesData
from .site import SiteData
from .plot import PlotData
from .tree import TreeData
from .model import ForestModel
from .io_utils import UVAFMEReader, UVAFMEWriter
from .climate import set_site_climate


class UVAFMEModel:
    """Main UVAFME model class."""
    
    def __init__(self):
        self.sites = []
        self.species_present = []
        self.filelist = ""
        self.reader = UVAFMEReader()
        self.writer = UVAFMEWriter()
        self.forest_model = ForestModel()
        
    def initialize_input_files(self, filelist=""):
        """Initialize input files and read data."""
        self.filelist = filelist
        
        # Read configuration files
        if filelist:
            file_config = self.reader.read_filelist(filelist)
            species_file = file_config.get('species_file', 'UVAFME2012_specieslist.csv')
            site_file = file_config.get('site_file', 'UVAFME2012_site.csv')
            climate_file = file_config.get('climate_file', 'UVAFME2012_climate.csv')
            config_file = file_config.get('config_file', 'uvafme_config.json')
        else:
            species_file = 'UVAFME2012_specieslist.csv'
            site_file = 'UVAFME2012_site.csv'
            climate_file = 'UVAFME2012_climate.csv'
            config_file = 'uvafme_config.json'
        
        # Load parameters
        try:
            params.load_from_file(config_file)
        except Exception as e:
            print(f"Using default parameters: {e}")
        
        # Read species data
        self.species_present = self.reader.read_species_file(species_file)
        print(f"Loaded {len(self.species_present)} species")
        
        # Read site data
        self.sites = self.reader.read_site_file(site_file)
        print(f"Loaded {len(self.sites)} sites")
        
        # Read climate data
        climate_data = self.reader.read_climate_file(climate_file)
        
        # Attach climate data to sites
        for site in self.sites:
            if site.site_id in climate_data:
                climate = climate_data[site.site_id]
                site.attach_climate(climate['tmin'], climate['tmax'], climate['precip'])
                site.attach_climate_std(climate['tmin_std'], climate['tmax_std'], climate['precip_std'])
            else:
                print(f"Warning: No climate data found for site {site.site_id}")
    
    def initialize_sitelist(self):
        """Initialize site list and species data."""
        for site in self.sites:
            # Attach species to site
            site.attach_species(self.species_present)
            
            # Initialize plots for each site
            site.numplots = params.numplots
            site.plots = []
            
            for plot_id in range(params.numplots):
                plot = PlotData()
                plot.initialize_plot(site.species, params.maxtrees, params.maxheight)
                site.plots.append(plot)
            
            # Initialize forest with starting trees
            self.forest_model.initialize_forest(site)
    
    def initialize_output_files(self, species_present):
        """Initialize output files."""
        self.writer.initialize_output_files(species_present)
    
    def draw_banner(self, numsites, species_present):
        """Draw the startup banner."""
        print("=" * 80)
        print("                       UVA Forest Model Enhanced")
        print("               Center For Regional Environmental Studies")
        print("                        University of Virginia")
        print("                   Department of Environmental Sciences")
        print("=" * 80)
        print()
        print("Running with parameters:")
        print(f"Number of sites: {numsites}")
        print(f"Number of years: {params.numyears}")
        print(f"Number of species: {len(species_present)}")
        print(f"Maximum number of trees: {params.maxtrees}")
        print(f"Maximum height of trees: {params.maxheight}")
        print(f"Plotsize: {params.plotsize}")
        print(f"Root depth: {params.rootdepth}")
        
        if params.with_clim_change:
            print("Running with climate change")
            print(f"Beginning at year: {params.begin_change_year}")
            print(f"Duration in years: {params.duration_of_change}")
        
        print(f"Printing interval in years: {params.year_print_interval}")
        print("=" * 80)
        print()
    
    def show_progress(self, site):
        """Show progress for current site."""
        num_site_species = len(site.species)
        print(f"Running for site {site.site_id} {site.site_name}")
        print(f"Number of species present {num_site_species}")
        
        if hasattr(site, 'altitude') and site.altitude != params.rnvalid:
            print(f"      Site altitude adjustment {site.altitude}")
        
        print(f"Site parameters: elevation {site.elevation:.2f}   slope     {site.slope:.2f}")
        print(f"                 fire/1000   {site.fire_prob:.2f}   wind/1000 {site.wind_prob:.2f}")
        print(f"                 SAFC*rtdpth {site.soil.A_field_cap:.2f}   A0_C      {site.soil.A0_c0:.2f}  A0_N {site.soil.A0_n0:.2f}")
        
    def bio_geo_climate(self, site, year):
        """Biogeochemical climate processing."""
        self.forest_model.bio_geo_climate(site, year)
    
    def write_soil_data(self, site, year):
        """Write soil data."""
        self.writer.write_soil_data(site, year)
    
    def write_site_data(self, site, year):
        """Write site data."""
        self.writer.write_site_data(site, year)
    
    def canopy(self, site):
        """Process canopy dynamics."""
        self.forest_model.canopy(site)
    
    def growth(self, site):
        """Process tree growth."""
        self.forest_model.growth(site)
    
    def mortality(self, site):
        """Process tree mortality."""
        self.forest_model.mortality(site)
    
    def renewal(self, site):
        """Process tree renewal."""
        self.forest_model.renewal(site)
    
    def write_genus_data(self, site, species_present, year):
        """Write genus data."""
        # Group species by genus
        genera = list(set(species.genus_name for species in species_present))
        
        for plot in site.plots:
            genus_data = plot.sum_over_sg(genera, field='genus')
            
            # Write genus-level data (implementation depends on output format)
            # This is a placeholder for genus-specific output
            pass
    
    def write_species_data(self, site, species_present, year):
        """Write species data."""
        self.writer.write_species_data(site, species_present, year)
    
    def write_tree_data(self, site, year):
        """Write tree data."""
        self.writer.write_tree_data(site, year)
    
    def close_output_files(self):
        """Close output files."""
        self.writer.close_output_files()
    
    def set_site_rng_seed(self, fixed_seed):
        """Set random number generator seed."""
        self.forest_model.set_site_rng_seed(fixed_seed)
    
    def run(self, filelist=""):
        """Main model execution loop."""
        
        # Handle command line arguments
        if len(sys.argv) > 1:
            filelist = sys.argv[1]
        
        # Prepare input files
        self.initialize_input_files(filelist)
        
        # Prepare site and species data
        self.initialize_sitelist()
        
        # Prepare output files
        self.initialize_output_files(self.species_present)
        
        # Write runtime vars to screen
        self.draw_banner(len(self.sites), self.species_present)
        
        # Start timing
        start_time = time.time()
        
        # Main site loop
        for sndx, site in enumerate(self.sites):
            
            # Check if site exists
            if site.site_wmo == params.rnvalid:
                print(f"             No site or climate file for site {site.site_id}")
                print(f"             Skipping site {site.site_name}")
                print()
                continue
            
            # Load climate and site specific vars, then adjust for altitude
            self.set_site_rng_seed(params.fixed_seed)
            set_site_climate(params.same_climate, params.fixed_seed)
            
            # Skip this site if no species are present
            if len(site.species) == 0:
                print(f"              No species present in site {site.site_id}")
                print(f"              Skipping site {site.site_name}")
                print()
                continue
            
            self.show_progress(site)
            
            # Run the model for this site
            for year in range(params.numyears + 1):
                
                # Run annual cycle
                self.forest_model.run_annual_cycle(site, year)
                
                # Write soil/CN/clim data for last year's trees but after BioGeo
                self.write_soil_data(site, year)
                self.write_site_data(site, year)
                
                # Determine print interval
                if params.spinup:
                    if year < params.spinup_yrs:
                        print_interval = 10 * params.year_print_interval
                    else:
                        print_interval = params.year_print_interval
                else:
                    print_interval = params.year_print_interval
                
                # Write output data
                if (year % print_interval == 0) or (year == params.numyears):
                    self.write_genus_data(site, self.species_present, year)
                    self.write_species_data(site, self.species_present, year)
                    if params.tree_level_data:
                        self.write_tree_data(site, year)
                
                # Progress indicator
                if year % 10 == 0:
                    stats = site.plots[0].get_statistics() if site.plots else {}
                    print(f"  Year {year}: {stats.get('live_trees', 0)} trees, "
                          f"BiomC={stats.get('total_biomass_c', 0):.2f}")
            
            # Report timing
            total_time = time.time()
            print(f"Cumulative time: {total_time - start_time:.2f} seconds")
            print("=" * 80)
        
        # Close output files
        self.close_output_files()


def main():
    """Main entry point for UVAFME model."""
    model = UVAFMEModel()
    model.run()


if __name__ == "__main__":
    main()