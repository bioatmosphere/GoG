"""
Parameter mapping between GAPPY and DEMENTpy models.

This module creates appropriate DEMENTpy runtime configurations and initialization
parameters based on GAPPY site and climate data.
"""

import numpy as np
from typing import Dict, Any, Optional


class ParameterMapper:
    """Maps GAPPY parameters to DEMENTpy initialization."""

    # Default DEMENTpy runtime parameters for forest ecosystems
    DEFAULT_DEMENT_CONFIG = {
        'pulse': 1,              # Number of pulses (years in this context)
        'end_time': 365,         # Days per pulse
        'interval': 30,          # Output interval (days)
        'dispersal': 0,          # Dispersal mode (0=default, 1=dispersal)
        'grid_size': 10,         # Grid dimensions (10x10 = 100 cells)
        'n_taxa': 5,             # Number of microbial taxa
        'n_enzymes': 3,          # Number of enzyme types
        'substrate_complexity': 3,  # Substrate complexity levels
    }

    # Forest soil initial conditions (typical temperate forest)
    DEFAULT_SOIL_INIT = {
        # Substrate pools (g/m²)
        'litter_C': 500.0,       # Initial litter C
        'litter_N': 20.0,        # Initial litter N (C:N = 25)
        'SOM_C': 5000.0,         # Soil organic matter C
        'SOM_N': 250.0,          # Soil organic matter N (C:N = 20)

        # Microbial biomass (g C/m²)
        'microbial_C': 50.0,     # Total microbial biomass
        'microbial_N': 5.0,      # Microbial N

        # Environmental conditions
        'moisture': 0.6,         # Relative moisture (0-1)
        'temperature': 15.0,     # °C
        'pH': 5.5,               # Soil pH
    }

    # Per-layer defaults for layered mode (AO organic, SA mineral-A, SB mineral-B)
    DEFAULT_LAYER_INIT = {
        'AO': {
            'substrate_C': 500.0,    # g/m² - litter layer
            'substrate_N': 16.7,     # g/m² (C:N ~ 30)
            'microbial_C': 15.0,     # 30% of total microbes
            'microbial_N': 1.67,
            'moisture': 0.5,         # Surface moisture
            'temp_dampen': 1.0,      # No attenuation at surface
        },
        'SA': {
            'substrate_C': 5000.0,   # g/m² - SOM layer
            'substrate_N': 250.0,    # g/m² (C:N ~ 20)
            'microbial_C': 25.0,     # 50% of total microbes
            'microbial_N': 2.78,
            'moisture': 0.7,         # Intermediate moisture
            'temp_dampen': 0.8,      # Moderate attenuation
        },
        'SB': {
            'substrate_C': 2000.0,   # g/m² - deep recalcitrant
            'substrate_N': 100.0,    # g/m² (C:N ~ 20)
            'microbial_C': 10.0,     # 20% of total microbes
            'microbial_N': 1.11,
            'moisture': 0.8,         # Deep moisture
            'temp_dampen': 0.6,      # Strong attenuation
        },
    }

    @staticmethod
    def create_dement_runtime(gappy_params: Optional[Dict] = None,
                             override: Optional[Dict] = None) -> Dict[str, Any]:
        """
        Create DEMENTpy runtime configuration from GAPPY parameters.

        Args:
            gappy_params: GAPPY model parameters (currently unused but for future)
            override: Dictionary of parameters to override defaults

        Returns:
            Complete DEMENTpy runtime configuration
        """
        config = ParameterMapper.DEFAULT_DEMENT_CONFIG.copy()

        # Apply GAPPY-specific mappings if provided
        if gappy_params:
            # Future: map GAPPY parameters to DEMENTpy
            # For example: adjust grid_size based on plot heterogeneity
            pass

        # Apply user overrides
        if override:
            config.update(override)

        return config

    @staticmethod
    def create_dement_initialization(gappy_site_data: Optional[Dict] = None,
                                    override: Optional[Dict] = None) -> Dict[str, Any]:
        """
        Create DEMENTpy initialization parameters from GAPPY site data.

        Args:
            gappy_site_data: Site data from GAPPY (soil C/N, climate)
            override: Dictionary to override defaults

        Returns:
            DEMENTpy initialization dictionary
        """
        init_params = ParameterMapper.DEFAULT_SOIL_INIT.copy()

        # Map GAPPY soil data if provided
        if gappy_site_data:
            # Map GAPPY soil C/N pools to DEMENTpy initial conditions
            if 'A0_c0' in gappy_site_data:
                # Convert from GAPPY units (tc/ha) to DEMENTpy units (g/m²)
                # 1 tc/ha = 100 g/m²
                init_params['litter_C'] = gappy_site_data['A0_c0'] * 100.0

            if 'A0_n0' in gappy_site_data:
                init_params['litter_N'] = gappy_site_data['A0_n0'] * 100.0

            if 'A_c0' in gappy_site_data:
                init_params['SOM_C'] = gappy_site_data['A_c0'] * 100.0

            if 'A_n0' in gappy_site_data:
                init_params['SOM_N'] = gappy_site_data['A_n0'] * 100.0

            # Map climate averages
            if 'mean_temp' in gappy_site_data:
                init_params['temperature'] = gappy_site_data['mean_temp']

            if 'mean_moisture' in gappy_site_data:
                init_params['moisture'] = gappy_site_data['mean_moisture']

        # Apply overrides
        if override:
            init_params.update(override)

        return init_params

    @staticmethod
    def create_layered_initialization(gappy_site_data: Optional[Dict] = None,
                                      override: Optional[Dict] = None) -> Dict[str, Dict[str, Any]]:
        """
        Create per-layer DEMENTpy initialization for the 3-grid layered mode.

        Maps GAPPY's A0_c0/n0, A_c0/n0, BL_c0/n0 to layer-specific init params.

        Args:
            gappy_site_data: Site data from GAPPY (soil C/N per layer)
            override: Per-layer overrides {layer_name: {param: value}}

        Returns:
            Dict with keys 'AO', 'SA', 'SB', each containing init params
        """
        import copy
        layers = copy.deepcopy(ParameterMapper.DEFAULT_LAYER_INIT)

        if gappy_site_data:
            # AO layer from A0 pools (1 tc/ha = 100 g/m²)
            if 'A0_c0' in gappy_site_data:
                layers['AO']['substrate_C'] = gappy_site_data['A0_c0'] * 100.0
            if 'A0_n0' in gappy_site_data:
                layers['AO']['substrate_N'] = gappy_site_data['A0_n0'] * 100.0

            # SA layer from A pools
            if 'A_c0' in gappy_site_data:
                layers['SA']['substrate_C'] = gappy_site_data['A_c0'] * 100.0
            if 'A_n0' in gappy_site_data:
                layers['SA']['substrate_N'] = gappy_site_data['A_n0'] * 100.0

            # SB layer from BL pools
            if 'BL_c0' in gappy_site_data:
                layers['SB']['substrate_C'] = gappy_site_data['BL_c0'] * 100.0
            if 'BL_n0' in gappy_site_data:
                layers['SB']['substrate_N'] = gappy_site_data['BL_n0'] * 100.0

            # Map climate if available
            if 'mean_temp' in gappy_site_data:
                for layer in layers.values():
                    layer['temperature'] = gappy_site_data['mean_temp']
            if 'mean_moisture' in gappy_site_data:
                layers['AO']['moisture'] = gappy_site_data['mean_moisture']
                layers['SA']['moisture'] = min(gappy_site_data['mean_moisture'] + 0.1, 1.0)
                layers['SB']['moisture'] = min(gappy_site_data['mean_moisture'] + 0.2, 1.0)

        if override:
            for layer_name, params in override.items():
                if layer_name in layers:
                    layers[layer_name].update(params)

        return layers

    @staticmethod
    def estimate_microbial_biomass(soil_c: float, soil_n: float) -> Dict[str, float]:
        """
        Estimate initial microbial biomass from soil C and N.

        Uses typical forest soil ratios:
        - Microbial C = 1-2% of total soil C
        - Microbial C:N ratio = 8:1 to 10:1

        Args:
            soil_c: Total soil C (g/m²)
            soil_n: Total soil N (g/m²)

        Returns:
            Dictionary with microbial_C and microbial_N
        """
        # Assume microbial C is 1.5% of total soil C
        microbial_c = soil_c * 0.015

        # Assume microbial C:N ratio of 9:1
        microbial_n = microbial_c / 9.0

        return {
            'microbial_C': microbial_c,
            'microbial_N': microbial_n,
            'microbial_CN': microbial_c / microbial_n
        }

    @staticmethod
    def create_spatial_config(n_plots: int, mode: str = 'aggregated') -> Dict[str, Any]:
        """
        Create spatial configuration for coupling.

        Args:
            n_plots: Number of GAPPY plots (typically 200)
            mode: 'aggregated', 'one_to_one', or 'layered'

        Returns:
            Spatial configuration dictionary
        """
        if mode == 'aggregated':
            return {
                'mode': 'aggregated',
                'n_dement_grids': 1,
                'aggregation_method': 'mean',  # How to aggregate litter inputs
                'plot_mapping': None  # No explicit mapping needed
            }
        elif mode == 'one_to_one':
            return {
                'mode': 'one_to_one',
                'n_dement_grids': n_plots,
                'aggregation_method': 'individual',
                'plot_mapping': list(range(n_plots))  # Direct 1:1 mapping
            }
        elif mode == 'layered':
            return {
                'mode': 'layered',
                'n_dement_grids': 3,
                'layer_names': ['AO', 'SA', 'SB'],
                'aggregation_method': 'cascading',
                'plot_mapping': None
            }
        else:
            raise ValueError(f"Unknown spatial mode: {mode}")

    @staticmethod
    def validate_parameters(config: Dict) -> bool:
        """
        Validate that DEMENTpy configuration is reasonable.

        Args:
            config: Configuration dictionary

        Returns:
            True if valid, raises ValueError otherwise
        """
        required_keys = ['pulse', 'end_time', 'interval']
        for key in required_keys:
            if key not in config:
                raise ValueError(f"Missing required parameter: {key}")

        if config['end_time'] <= 0:
            raise ValueError("end_time must be positive")

        if config['interval'] <= 0 or config['interval'] > config['end_time']:
            raise ValueError(f"interval must be between 1 and {config['end_time']}")

        return True

    @staticmethod
    def print_parameter_summary(runtime: Dict, init: Dict, spatial: Dict):
        """Print a summary of coupling parameters."""
        print("\n" + "="*60)
        print("DEMENTpy Coupling Configuration Summary")
        print("="*60)

        print("\nRuntime Parameters:")
        for key, value in runtime.items():
            print(f"  {key:20s}: {value}")

        print("\nInitialization Parameters:")
        for key, value in init.items():
            if isinstance(value, float):
                print(f"  {key:20s}: {value:.2f}")
            else:
                print(f"  {key:20s}: {value}")

        print("\nSpatial Configuration:")
        for key, value in spatial.items():
            if key != 'plot_mapping':  # Don't print long arrays
                print(f"  {key:20s}: {value}")

        print("="*60 + "\n")


def create_default_coupling_config(n_plots: int = 200,
                                  spatial_mode: str = 'aggregated') -> Dict[str, Any]:
    """
    Create a complete default coupling configuration.

    Args:
        n_plots: Number of GAPPY plots
        spatial_mode: 'aggregated', 'one_to_one', or 'layered'

    Returns:
        Complete configuration dictionary with runtime, init, and spatial settings
    """
    mapper = ParameterMapper()

    config = {
        'runtime': mapper.create_dement_runtime(),
        'initialization': mapper.create_dement_initialization(),
        'spatial': mapper.create_spatial_config(n_plots, spatial_mode)
    }

    # Add per-layer initialization for layered mode
    if spatial_mode == 'layered':
        config['layer_initialization'] = mapper.create_layered_initialization()

    # Validate
    mapper.validate_parameters(config['runtime'])

    return config


if __name__ == "__main__":
    # Demonstrate parameter mapping
    print("Creating default coupling configuration...\n")

    # Aggregated mode
    config_agg = create_default_coupling_config(n_plots=200, spatial_mode='aggregated')
    mapper = ParameterMapper()
    mapper.print_parameter_summary(
        config_agg['runtime'],
        config_agg['initialization'],
        config_agg['spatial']
    )

    # Example with GAPPY data
    print("\nExample with GAPPY site data:")
    gappy_site = {
        'A0_c0': 5.0,     # tc/ha (surface layer)
        'A0_n0': 0.2,     # tn/ha
        'A_c0': 50.0,     # tc/ha (mineral layer)
        'A_n0': 2.5,      # tn/ha
        'mean_temp': 12.5,  # °C
        'mean_moisture': 0.65
    }

    init_params = mapper.create_dement_initialization(gappy_site)
    print(f"\nMapped initialization parameters:")
    for key, val in init_params.items():
        if isinstance(val, float):
            print(f"  {key:20s}: {val:.2f}")
