"""
Model module for UVAFME vegetation model.
Orchestrates biogeochemical processes, forest dynamics, and ecosystem interactions.
"""

import random
import math
import numpy as np
from typing import List, Dict, Optional
from .constants import *
from .parameters import params
from .species import SpeciesData
from .site import SiteData
from .tree import TreeData
from .plot import PlotData
from .climate import ex_rad, hargrea, cov365, cov365a
from .soil import SoilData


class ForestModel:
    """Main forest dynamics model orchestrating all processes."""
    
    def __init__(self):
        self.random_seed = None
        self.climate_uniform_func = random.random
    
    def set_site_rng_seed(self, fixed_seed: bool = False):
        """Set random number generator seed for site."""
        if fixed_seed:
            self.random_seed = 42
            random.seed(self.random_seed)
        else:
            self.random_seed = None
            random.seed()
    
    def bio_geo_climate(self, site: SiteData, year: int):
        """Process biogeochemical climate interactions."""
        # Calculate daily climate variables
        self.calculate_daily_climate(site, year)
        
        # Calculate growing degree days
        self.calculate_degree_days(site)
        
        # Calculate potential evapotranspiration
        self.calculate_potential_evapotranspiration(site)
        
        # Update species climate response factors
        self.update_species_climate_responses(site)
        
        # Process soil biogeochemistry
        self.process_soil_biogeochemistry(site)
    
    def calculate_daily_climate(self, site: SiteData, year: int):
        """Calculate daily climate variables from monthly data."""
        # Convert monthly temperature to daily
        daily_tmin = cov365(site.tmin)
        daily_tmax = cov365(site.tmax)
        
        # Convert monthly precipitation to daily (with randomness)
        daily_precip = cov365a(site.precip, self.climate_uniform_func)
        
        # Calculate daily averages for the year
        site.rain = np.mean(daily_precip) / 10.0  # Convert to cm/day
        
        # Calculate temperature statistics
        daily_temp_avg = (daily_tmin + daily_tmax) / 2.0
        site.freeze = np.sum(daily_temp_avg < 0) / len(daily_temp_avg)
    
    def calculate_degree_days(self, site: SiteData):
        """Calculate growing degree days."""
        # Base temperature for growing degree days
        base_temp = 5.0
        
        # Calculate from monthly data
        monthly_avg = (site.tmin + site.tmax) / 2.0
        
        # Sum degree days above base temperature
        degree_days = np.sum(np.maximum(monthly_avg - base_temp, 0) * 30)  # 30 days/month
        site.deg_days = degree_days
        
        # Calculate growing days (days above base temperature)
        site.grow_days = np.sum(monthly_avg > base_temp) * 30
    
    def calculate_potential_evapotranspiration(self, site: SiteData):
        """Calculate potential evapotranspiration using Hargreaves method."""
        # Use middle of year for radiation calculation
        julian_day = 180  # Mid-year approximation
        
        # Calculate extraterrestrial radiation
        erad, daylength, exradmx = ex_rad(julian_day, site.latitude)
        
        # Calculate monthly potential evapotranspiration
        monthly_pet = []
        for month in range(12):
            tmin = site.tmin[month]
            tmax = site.tmax[month]
            tavg = (tmin + tmax) / 2.0
            
            # Hargreaves equation
            pet = hargrea(tmin, tmax, tavg, erad)
            monthly_pet.append(pet)
        
        # Average daily potential evapotranspiration
        site.pot_evap_day = np.mean(monthly_pet)
    
    def update_species_climate_responses(self, site: SiteData):
        """Update species climate response factors."""
        for plot in site.plots:
            for species in plot.species:
                # Temperature response
                species.temp_rsp(site.deg_days)
                
                # Drought response (simplified)
                dry_days = max(0, 100 - site.rain * 365)  # Estimate dry days
                dry_days_surface = dry_days * 0.8  # Surface drying faster
                species.drought_rsp(dry_days, dry_days_surface)
                
                # Flood response (simplified)
                flood_days = max(0, site.rain * 365 - 1000) / 100  # Estimate flood days
                species.flood_rsp(flood_days)
                
                # Fire response (based on site probability)
                fire_occurrence = 1 if random.random() < site.fire_prob else 0
                species.fire_rsp(fire_occurrence)
    
    def process_soil_biogeochemistry(self, site: SiteData):
        """Process soil biogeochemical cycles."""
        # Calculate litter inputs from existing trees
        total_litter_c1 = 0.0  # Above-ground litter
        total_litter_c2 = 0.0  # Below-ground litter
        total_litter_n1 = 0.0  # Above-ground N
        total_litter_n2 = 0.0  # Below-ground N
        
        for plot in site.plots:
            for tree in plot.trees:
                if tree.mort_marker:
                    # Dead tree contributes to litter
                    total_litter_c1 += tree.biomC * 0.7  # 70% above-ground
                    total_litter_c2 += tree.biomC * 0.3  # 30% below-ground
                    total_litter_n1 += tree.biomN * 0.7
                    total_litter_n2 += tree.biomN * 0.3
                else:
                    # Living tree contributes small amount (leaf fall, etc.)
                    annual_turnover = 0.1  # 10% annual turnover
                    total_litter_c1 += tree.leaf_bm * annual_turnover
                    total_litter_n1 += tree.leaf_bm * annual_turnover / DEC_LEAF_C_N
        
        # Average temperature for soil processes
        temp_avg = np.mean((site.tmin + site.tmax) / 2.0)
        
        # Soil water factors (simplified)
        aow0_scaled = 0.5  # Placeholder
        saw0_scaled = 0.7  # Placeholder
        sbw0_scaled = 0.8  # Placeholder
        
        # Process soil decomposition
        avail_N, C_resp = site.soil.soil_decomp(
            total_litter_c1, total_litter_c2, total_litter_n1, total_litter_n2,
            temp_avg, site.rain, aow0_scaled, saw0_scaled, sbw0_scaled
        )
        
        # Process soil water balance
        water_results = site.soil.soil_water(
            site.slope, site.leaf_area_ind, site.leaf_area_w0,
            site.sigma, site.freeze, site.rain, site.pot_evap_day
        )
        
        site.act_evap_day = water_results[0]
        
        # Calculate dry days from soil water
        site.dry_days_upper_layer = max(0, 100 - water_results[3] * 100)
        site.dry_days_base_layer = max(0, 100 - water_results[5] * 100)
        
        # Calculate flood days from soil water
        site.flood_days = max(0, water_results[1] * 100 - 80)
    
    def canopy(self, site: SiteData):
        """Process canopy light interactions."""
        for plot in site.plots:
            # Calculate light profile through canopy
            plot.calculate_light_profile()
            
            # Update tree light responses
            for tree in plot.trees:
                # Get light availability at tree height (default to moderate light)
                height_idx = min(int(tree.forska_ht), len(plot.con_light) - 1)
                
                if tree.conifer:
                    available_light = plot.con_light[height_idx] if height_idx < len(plot.con_light) and height_idx >= 0 else 0.6
                else:
                    available_light = plot.dec_light[height_idx] if height_idx < len(plot.dec_light) and height_idx >= 0 else 0.6
                
                # Apply light response to tree
                light_factor = tree.light_rsp(available_light)
                
                # Update tree's environmental stress
                tree.env_stress_factor = tree.env_stress(light_factor)
    
    def growth(self, site: SiteData):
        """Process tree growth."""
        for plot in site.plots:
            for tree in plot.trees:
                if not tree.mort_marker:
                    # Calculate maximum growth potential
                    tree.max_growth()
                    
                    # Apply environmental stress to growth
                    env_stress = getattr(tree, 'env_stress_factor', 0.5)  # Default moderate stress
                    actual_growth = tree.diam_max * env_stress
                    
                    # Update tree diameter
                    tree.diam_bht += actual_growth
                    
                    # Ensure diameter doesn't exceed maximum
                    tree.diam_bht = min(tree.diam_bht, tree.max_diam)
                    
                    # Update tree metrics
                    tree.calculate_all_metrics()
                    
                    # Check for growth stress
                    if actual_growth < tree.diam_max * 0.1:  # Less than 10% of potential
                        tree.mort_marker = True
    
    def mortality(self, site: SiteData):
        """Process tree mortality."""
        for plot in site.plots:
            for tree in plot.trees:
                if not tree.mort_marker:
                    # Check age-based mortality
                    if not tree.age_survival():
                        tree.mort_marker = True
                        continue
                    
                    # Check growth-based mortality
                    if not tree.growth_survival():
                        tree.mort_marker = True
                        continue
                    
                    # Environmental mortality (fire, wind, etc.)
                    if random.random() < site.fire_prob * (1.0 - tree.fc_fire):
                        tree.mort_marker = True
                        continue
                    
                    if random.random() < site.wind_prob * (1.0 - tree.fc_wind):
                        tree.mort_marker = True
                        continue
    
    def renewal(self, site: SiteData):
        """Process tree renewal (seedling establishment and growth)."""
        for plot in site.plots:
            # Remove dead trees
            plot.remove_dead_trees()
            
            # Calculate seedling establishment
            for i, species in enumerate(plot.species):
                # Calculate seed production from existing trees
                seed_production = 0.0
                for tree in plot.trees:
                    if tree.unique_id == species.unique_id and not tree.mort_marker:
                        # Seed production based on tree size and species parameters
                        seed_production += tree.diam_bht * species.seed_num / 100.0
                
                # Calculate seedling establishment
                establishment_rate = seed_production * species.seed_surv * species.seedling_lg
                
                # Environmental filters
                establishment_rate *= species.fc_degday * species.fc_drought * species.fc_flood
                
                # Add to seedling bank
                plot.seedling[i] += establishment_rate
                
                # Promote seedlings to saplings
                new_saplings = int(plot.seedling[i] * 0.1)  # 10% promotion rate
                plot.seedling[i] -= new_saplings
                
                # Create new trees from saplings
                for _ in range(new_saplings):
                    if len(plot.trees) < params.maxtrees:
                        new_tree = TreeData()
                        new_tree.initialize_tree(species, i)
                        
                        # Initialize small tree
                        new_tree.diam_bht = 1.0 + random.random() * 2.0  # 1-3 cm
                        new_tree.calculate_all_metrics()
                        
                        plot.add_tree(new_tree)
    
    def initialize_forest(self, site: SiteData):
        """Initialize forest with starting trees."""
        for plot in site.plots:
            # Add some initial trees of different species
            for i, species in enumerate(plot.species):
                # Number of initial trees (minimum 5 per species)
                n_initial = max(5, int(10 * species.invader + 10))
                
                for _ in range(min(n_initial, params.maxtrees // len(plot.species))):
                    tree = TreeData()
                    tree.initialize_tree(species, i)
                    
                    # Initialize with random sizes
                    tree.diam_bht = random.uniform(5.0, 30.0)  # 5-30 cm DBH
                    tree.calculate_all_metrics()
                    
                    plot.add_tree(tree)
                    
            print(f"  Initialized plot with {len(plot.trees)} trees")
    
    def run_annual_cycle(self, site: SiteData, year: int):
        """Run complete annual cycle for a site."""
        # Biogeochemical processes
        self.bio_geo_climate(site, year)
        
        # Forest dynamics
        self.canopy(site)
        self.growth(site)
        self.mortality(site)
        self.renewal(site)
        
        # Update site-level statistics
        self.update_site_statistics(site)
    
    def update_site_statistics(self, site: SiteData):
        """Update site-level forest statistics."""
        total_biomass_c = 0.0
        total_biomass_n = 0.0
        total_basal_area = 0.0
        total_trees = 0
        
        for plot in site.plots:
            plot_stats = plot.get_statistics()
            total_biomass_c += plot_stats['total_biomass_c']
            total_biomass_n += plot_stats['total_biomass_n']
            total_basal_area += plot_stats['total_basal_area']
            total_trees += plot_stats['total_trees']
        
        # Store site-level statistics
        site.total_biomass_c = total_biomass_c
        site.total_biomass_n = total_biomass_n
        site.total_basal_area = total_basal_area
        site.total_trees = total_trees