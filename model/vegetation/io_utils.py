"""
Input/Output utilities for UVAFME vegetation model.
Handles reading and writing of climate, species, and site data.
"""

import os
import csv
import json
import numpy as np
from typing import List, Dict, Tuple, Optional
from .constants import *
from .species import SpeciesData
from .site import SiteData
from .parameters import Parameters


class UVAFMEReader:
    """Handles reading of UVAFME input files."""
    
    def __init__(self, base_path: str = "input_data"):
        self.base_path = base_path
        self.ensure_directories()
    
    def ensure_directories(self):
        """Ensure input directories exist."""
        os.makedirs(self.base_path, exist_ok=True)
    
    def read_species_file(self, filename: str = "species.csv") -> List[SpeciesData]:
        """Read species data from CSV file."""
        filepath = os.path.join(self.base_path, filename)
        species_list = []
        
        if not os.path.exists(filepath):
            print(f"Warning: Species file {filepath} not found, creating example")
            self.create_example_species_file(filepath)
        
        with open(filepath, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                species = SpeciesData()
                species.initialize_species(
                    species_id=int(row['species_id']),
                    genus_name=row['genus_name'],
                    taxonomic_name=row['taxonomic_name'],
                    unique_id=row['unique_id'],
                    common_name=row['common_name'],
                    genus_id=int(row['genus_id']),
                    shade_tol=int(row['shade_tol']),
                    lownutr_tol=int(row['lownutr_tol']),
                    stress_tol=int(row['stress_tol']),
                    age_tol=int(row['age_tol']),
                    drought_tol=int(row['drought_tol']),
                    flood_tol=int(row['flood_tol']),
                    fire_tol=int(row['fire_tol']),
                    max_age=float(row['max_age']),
                    max_diam=float(row['max_diam']),
                    max_ht=float(row['max_ht']),
                    wood_bulk_dens=float(row['wood_bulk_dens']),
                    rootdepth=float(row['rootdepth']),
                    leafdiam_a=float(row['leafdiam_a']),
                    leafarea_c=float(row['leafarea_c']),
                    deg_day_min=float(row['deg_day_min']),
                    deg_day_opt=float(row['deg_day_opt']),
                    deg_day_max=float(row['deg_day_max']),
                    seedling_lg=float(row['seedling_lg']),
                    invader=float(row['invader']),
                    seed_num=float(row['seed_num']),
                    sprout_num=float(row['sprout_num']),
                    seed_surv=float(row['seed_surv']),
                    arfa_0=float(row['arfa_0']),
                    g=float(row['g']),
                    conifer=row['conifer'].lower() == 'true'
                )
                species_list.append(species)
        
        return species_list
    
    def read_site_file(self, filename: str = "sites.csv") -> List[SiteData]:
        """Read site data from CSV file."""
        filepath = os.path.join(self.base_path, filename)
        sites_list = []
        
        if not os.path.exists(filepath):
            print(f"Warning: Site file {filepath} not found, creating example")
            self.create_example_site_file(filepath)
        
        with open(filepath, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                site = SiteData()
                
                # Parse temperature and precipitation lapse rates
                temp_lapse = [float(x) for x in row['temp_lapse'].split(',')]
                prcp_lapse = [float(x) for x in row['prcp_lapse'].split(',')]
                
                site.initialize_site(
                    siteid=int(row['site_id']),
                    sitename=row['site_name'],
                    siteregion=row['region'],
                    lat=float(row['latitude']),
                    long=float(row['longitude']),
                    wmo=float(row['site_wmo']),
                    elevation=float(row['elevation']),
                    slope=float(row['slope']),
                    Afc=float(row['A_field_cap']),
                    A_perm_wp=float(row['A_perm_wp']),
                    lai=float(row['leaf_area_ind']),
                    base_h=float(row['base_h']),
                    lai_w0=float(row['leaf_area_w0']),
                    A0_w0=float(row['A0_w0']),
                    A_w0=float(row['A_w0']),
                    sbase_w0=float(row['sbase_w0']),
                    fire_prob=float(row['fire_prob']),
                    wind_prob=float(row['wind_prob']),
                    A0_c0=float(row['A0_c0']),
                    A0_n0=float(row['A0_n0']),
                    A_c0=float(row['A_c0']),
                    A_n0=float(row['A_n0']),
                    sbase_c0=float(row['sbase_c0']),
                    sbase_n0=float(row['sbase_n0']),
                    sigma=float(row['sigma']),
                    temp_lapse=temp_lapse,
                    prcp_lapse=prcp_lapse
                )
                sites_list.append(site)
        
        return sites_list
    
    def read_climate_file(self, filename: str = "climate.csv") -> Dict[int, Dict[str, np.ndarray]]:
        """Read climate data from CSV file."""
        filepath = os.path.join(self.base_path, filename)
        climate_data = {}
        
        if not os.path.exists(filepath):
            print(f"Warning: Climate file {filepath} not found, creating example")
            self.create_example_climate_file(filepath)
        
        with open(filepath, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                site_id = int(row['site_id'])
                
                # Parse monthly data
                tmin = [float(row[f'tmin_{i}']) for i in range(1, 13)]
                tmax = [float(row[f'tmax_{i}']) for i in range(1, 13)]
                precip = [float(row[f'precip_{i}']) for i in range(1, 13)]
                
                # Optional standard deviations
                tmin_std = [float(row.get(f'tmin_std_{i}', 0)) for i in range(1, 13)]
                tmax_std = [float(row.get(f'tmax_std_{i}', 0)) for i in range(1, 13)]
                precip_std = [float(row.get(f'precip_std_{i}', 0)) for i in range(1, 13)]
                
                climate_data[site_id] = {
                    'tmin': np.array(tmin),
                    'tmax': np.array(tmax),
                    'precip': np.array(precip),
                    'tmin_std': np.array(tmin_std),
                    'tmax_std': np.array(tmax_std),
                    'precip_std': np.array(precip_std)
                }
        
        return climate_data
    
    def read_filelist(self, filename: str = "filelist.txt") -> Dict[str, str]:
        """Read file list configuration."""
        filepath = os.path.join(self.base_path, filename)
        filelist = {}
        
        if not os.path.exists(filepath):
            print(f"Warning: Filelist {filepath} not found, creating example")
            self.create_example_filelist(filepath)
        
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    key, value = line.split('=', 1)
                    filelist[key.strip()] = value.strip()
        
        return filelist
    
    def create_example_species_file(self, filepath: str):
        """Create example species file."""
        example_species = [
            {
                'species_id': 1, 'genus_name': 'Quercus', 'taxonomic_name': 'Quercus alba',
                'unique_id': 'QUAL', 'common_name': 'White Oak', 'genus_id': 1,
                'shade_tol': 3, 'lownutr_tol': 2, 'stress_tol': 2, 'age_tol': 5,
                'drought_tol': 2, 'flood_tol': 2, 'fire_tol': 2, 'max_age': 300.0,
                'max_diam': 150.0, 'max_ht': 35.0, 'wood_bulk_dens': 0.68,
                'rootdepth': 2.0, 'leafdiam_a': 0.1, 'leafarea_c': 50.0,
                'deg_day_min': 800.0, 'deg_day_opt': 1200.0, 'deg_day_max': 1800.0,
                'seedling_lg': 0.5, 'invader': 0.0, 'seed_num': 100.0,
                'sprout_num': 50.0, 'seed_surv': 0.1, 'arfa_0': 0.8, 'g': 0.3,
                'conifer': 'false'
            },
            {
                'species_id': 2, 'genus_name': 'Pinus', 'taxonomic_name': 'Pinus strobus',
                'unique_id': 'PIST', 'common_name': 'White Pine', 'genus_id': 2,
                'shade_tol': 2, 'lownutr_tol': 3, 'stress_tol': 3, 'age_tol': 4,
                'drought_tol': 3, 'flood_tol': 1, 'fire_tol': 1, 'max_age': 250.0,
                'max_diam': 120.0, 'max_ht': 40.0, 'wood_bulk_dens': 0.35,
                'rootdepth': 1.5, 'leafdiam_a': 0.15, 'leafarea_c': 40.0,
                'deg_day_min': 600.0, 'deg_day_opt': 1000.0, 'deg_day_max': 1600.0,
                'seedling_lg': 0.3, 'invader': 0.0, 'seed_num': 200.0,
                'sprout_num': 0.0, 'seed_surv': 0.05, 'arfa_0': 0.9, 'g': 0.4,
                'conifer': 'true'
            }
        ]
        
        with open(filepath, 'w', newline='') as csvfile:
            fieldnames = list(example_species[0].keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(example_species)
    
    def create_example_site_file(self, filepath: str):
        """Create example site file."""
        example_site = {
            'site_id': 1, 'site_name': 'Test Site', 'region': 'Virginia',
            'latitude': 37.0, 'longitude': -78.0, 'site_wmo': 123456,
            'elevation': 200.0, 'slope': 5.0, 'A_field_cap': 20.0,
            'A_perm_wp': 8.0, 'leaf_area_ind': 3.0, 'base_h': 70.0,
            'leaf_area_w0': 0.5, 'A0_w0': 2.0, 'A_w0': 15.0, 'sbase_w0': 40.0,
            'fire_prob': 0.01, 'wind_prob': 0.05, 'A0_c0': 10.0, 'A0_n0': 1.0,
            'A_c0': 50.0, 'A_n0': 5.0, 'sbase_c0': 100.0, 'sbase_n0': 10.0,
            'sigma': 0.1, 'temp_lapse': ','.join(['0.6'] * 12),
            'prcp_lapse': ','.join(['0.1'] * 12)
        }
        
        with open(filepath, 'w', newline='') as csvfile:
            fieldnames = list(example_site.keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(example_site)
    
    def create_example_climate_file(self, filepath: str):
        """Create example climate file."""
        # Example climate data for Virginia
        tmin = [-5, -2, 3, 8, 13, 18, 20, 19, 15, 9, 4, -1]
        tmax = [5, 8, 15, 22, 27, 32, 35, 33, 28, 22, 14, 7]
        precip = [80, 70, 90, 100, 120, 110, 90, 85, 75, 65, 70, 85]
        
        example_climate = {'site_id': 1}
        
        for i, (tmin_val, tmax_val, precip_val) in enumerate(zip(tmin, tmax, precip), 1):
            example_climate[f'tmin_{i}'] = tmin_val
            example_climate[f'tmax_{i}'] = tmax_val
            example_climate[f'precip_{i}'] = precip_val
        
        with open(filepath, 'w', newline='') as csvfile:
            fieldnames = list(example_climate.keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(example_climate)
    
    def create_example_filelist(self, filepath: str):
        """Create example filelist file."""
        example_filelist = """# UVAFME Input File List
# Format: key=value
species_file=species.csv
site_file=sites.csv
climate_file=climate.csv
config_file=uvafme_config.json
"""
        with open(filepath, 'w') as f:
            f.write(example_filelist)


class UVAFMEWriter:
    """Handles writing of UVAFME output files."""
    
    def __init__(self, base_path: str = "output_data"):
        self.base_path = base_path
        self.ensure_directories()
        self.output_files = {}
    
    def ensure_directories(self):
        """Ensure output directories exist."""
        os.makedirs(self.base_path, exist_ok=True)
    
    def initialize_output_files(self, species_present=None):
        """Initialize output files for writing."""
        # Species-level output
        species_file = os.path.join(self.base_path, "species_data.csv")
        self.output_files['species'] = open(species_file, 'w', newline='')
        species_writer = csv.writer(self.output_files['species'])
        species_writer.writerow(['year', 'site_id', 'species_id', 'genus_name', 
                                'biomass_c', 'biomass_n', 'basal_area', 'n_trees'])
        
        # Site-level output
        site_file = os.path.join(self.base_path, "site_data.csv")
        self.output_files['site'] = open(site_file, 'w', newline='')
        site_writer = csv.writer(self.output_files['site'])
        site_writer.writerow(['year', 'site_id', 'rain', 'pot_evap', 'act_evap',
                             'grow_days', 'deg_days', 'dry_days_upper', 'dry_days_base',
                             'flood_days'])
        
        # Soil-level output
        soil_file = os.path.join(self.base_path, "soil_data.csv")
        self.output_files['soil'] = open(soil_file, 'w', newline='')
        soil_writer = csv.writer(self.output_files['soil'])
        soil_writer.writerow(['year', 'site_id', 'A0_c0', 'A_c0', 'A0_n0', 'A_n0',
                             'BL_c0', 'BL_n0', 'total_C_rsp', 'avail_N'])
    
    def write_species_data(self, site, species_present, year):
        """Write species-level data."""
        if 'species' not in self.output_files:
            return
        
        writer = csv.writer(self.output_files['species'])
        for plot in site.plots:
            for species in plot.species:
                # Get species metrics from plot
                metrics = plot.sum_over_sg([species.unique_id], field='species')
                
                writer.writerow([
                    year, site.site_id, species.species_id, species.genus_name,
                    metrics['biomC'][0], metrics['biomN'][0], 
                    metrics['basal_area'][0], metrics['n_trees'][0]
                ])
    
    def write_site_data(self, site, year):
        """Write site-level data."""
        if 'site' not in self.output_files:
            return
        
        writer = csv.writer(self.output_files['site'])
        writer.writerow([
            year, site.site_id, site.rain, site.pot_evap_day, site.act_evap_day,
            site.grow_days, site.deg_days, site.dry_days_upper_layer,
            site.dry_days_base_layer, site.flood_days
        ])
    
    def write_soil_data(self, site, year):
        """Write soil-level data."""
        if 'soil' not in self.output_files:
            return
        
        writer = csv.writer(self.output_files['soil'])
        soil = site.soil
        writer.writerow([
            year, site.site_id, soil.A0_c0, soil.A_c0, soil.A0_n0, soil.A_n0,
            soil.BL_c0, soil.BL_n0, soil.total_C_rsp, soil.avail_N
        ])
    
    def write_tree_data(self, site, year):
        """Write tree-level data (if enabled)."""
        tree_file = os.path.join(self.base_path, f"trees_year_{year}.csv")
        
        with open(tree_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['year', 'site_id', 'plot_id', 'tree_id', 'species_id',
                            'genus_name', 'dbh', 'height', 'biomass_c', 'biomass_n',
                            'leaf_biomass', 'mort_marker'])
            
            for plot_id, plot in enumerate(site.plots):
                for tree_id, tree in enumerate(plot.trees):
                    writer.writerow([
                        year, site.site_id, plot_id, tree_id, tree.species_id,
                        tree.genus_name, tree.diam_bht, tree.forska_ht,
                        tree.biomC, tree.biomN, tree.leaf_bm, tree.mort_marker
                    ])
    
    def close_output_files(self):
        """Close all output files."""
        for file_handle in self.output_files.values():
            file_handle.close()
        self.output_files.clear()