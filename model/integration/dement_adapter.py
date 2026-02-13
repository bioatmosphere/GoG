"""
DEMENTpy Adapter - Coupling layer between GAPPY and DEMENTpy.

This module provides the main adapter class that replaces GAPPY's empirical
soil_decomp() function with mechanistic DEMENTpy microbial decomposition.
"""

import sys
import os
import numpy as np
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path

# Handle both package and script imports
try:
    from .unit_conversions import UnitConverter
    from .parameter_mapping import ParameterMapper, create_default_coupling_config
except ImportError:
    from unit_conversions import UnitConverter
    from parameter_mapping import ParameterMapper, create_default_coupling_config

# Import DEMENTpy library interface
DEMENT_AVAILABLE = False
DEMENTLibrary = None

try:
    # Add DEMENTpy src to path
    dementpy_src = Path(__file__).parent.parent / 'microbiome' / 'DEMENTpy' / 'src'
    if str(dementpy_src) not in sys.path:
        sys.path.insert(0, str(dementpy_src))

    # Try to import - this will test if pandas and other deps are available
    from library_interface import DEMENTLibrary
    DEMENT_AVAILABLE = True
    print(f"✓ DEMENTpy library_interface imported successfully")
    print(f"  Source: {dementpy_src}")

except ImportError as e:
    print(f"ℹ️  DEMENTpy library_interface not available: {e}")
    print(f"  Adapter will use fallback empirical model")
    print(f"  (This is normal if pandas is not in the current environment)")
except Exception as e:
    print(f"⚠️  Unexpected error importing DEMENTpy: {e}")
    print(f"  Adapter will use fallback empirical model")


class DEMENTpyAdapter:
    """
    Adapter for coupling GAPPY vegetation model with DEMENTpy microbial model.

    This class manages one or more DEMENTpy Grid instances and provides an
    interface that matches GAPPY's soil_decomp() function signature.

    Attributes:
        spatial_mode: Either 'aggregated' (1 grid) or 'one_to_one' (200 grids)
        n_grids: Number of DEMENTpy grid instances
        dement_grids: List of DEMENTpy Grid instances (placeholder for now)
        runtime_config: DEMENTpy runtime configuration
        converter: Unit conversion utility
        current_year: Current simulation year
        state_initialized: Whether DEMENTpy state is initialized
    """

    def __init__(self,
                 n_plots: int = 200,
                 spatial_mode: str = 'aggregated',
                 runtime_config: Optional[Dict] = None,
                 init_config: Optional[Dict] = None,
                 enable_dement: bool = False):
        """
        Initialize the DEMENTpy adapter.

        Args:
            n_plots: Number of GAPPY plots (typically 200)
            spatial_mode: 'aggregated', 'one_to_one', or 'layered'
            runtime_config: DEMENTpy runtime parameters (optional)
            init_config: DEMENTpy initialization parameters (optional)
            enable_dement: If True, actually use DEMENTpy; if False, use fallback
        """
        self.n_plots = n_plots
        self.spatial_mode = spatial_mode

        # Check if DEMENTpy is available before enabling
        if enable_dement and not DEMENT_AVAILABLE:
            print("Warning: DEMENTpy requested but not available. Using fallback model.")
            enable_dement = False

        self.enable_dement = enable_dement
        self.converter = UnitConverter()
        self.mapper = ParameterMapper()

        # Create configuration
        full_config = create_default_coupling_config(n_plots, spatial_mode)

        if runtime_config:
            full_config['runtime'].update(runtime_config)
        if init_config:
            full_config['initialization'].update(init_config)

        self.runtime_config = full_config['runtime']
        self.init_config = full_config['initialization']
        self.spatial_config = full_config['spatial']

        # Layer initialization for layered mode
        self.layer_init = full_config.get('layer_initialization', None)

        # Determine number of grids to create
        if spatial_mode == 'aggregated':
            self.n_grids = 1
        elif spatial_mode == 'one_to_one':
            self.n_grids = n_plots
        elif spatial_mode == 'layered':
            self.n_grids = 3
        else:
            raise ValueError(f"Unknown spatial_mode: {spatial_mode}")

        # Layer constants for empirical transfer rates (layered mode)
        self.AO_CN_0 = 30.0
        self.SA_CN_0 = 4.0
        self.SB_CN_0 = 20.0
        self.AO_RESP = 5.24e-4
        self.SA_RESP = 1.24e-5
        self.SB_RESP = 2.74e-7

        # Per-layer C/N pools in GAPPY units (tc/ha, tn/ha) — layered mode only
        self.layer_pools = {
            'ao_c': 0.0, 'ao_n': 0.0,
            'sa_c': 0.0, 'sa_n': 0.0,
            'sb_c': 0.0, 'sb_n': 0.0,
        }

        # Initialize DEMENTpy grids (placeholder for now)
        self.dement_grids = []
        self.state_initialized = False
        self.current_year = 0
        self.current_day = 0  # Track current day within year (0-364)

        # Statistics tracking
        self.coupling_stats = {
            'total_calls': 0,
            'total_litter_c': 0.0,
            'total_litter_n': 0.0,
            'total_resp': 0.0,
            'total_n_avail': 0.0,
        }

        print(f"DEMENTpyAdapter initialized:")
        print(f"  Spatial mode: {spatial_mode}")
        print(f"  Number of grids: {self.n_grids}")
        print(f"  DEMENTpy enabled: {enable_dement}")
        if not enable_dement:
            print(f"  -> Using fallback empirical model")

    def initialize_dement_grids(self, gappy_site_data: Optional[Dict] = None):
        """
        Initialize DEMENTpy Grid instances.

        Creates actual DEMENTpy Grid instances using the library interface.
        In layered mode, creates 3 grids with layer-specific initial conditions.

        Args:
            gappy_site_data: Site data from GAPPY for initialization
        """
        if self.state_initialized:
            print("DEMENTpy grids already initialized")
            return

        # Update initialization config with GAPPY data if provided
        if gappy_site_data:
            init_params = self.mapper.create_dement_initialization(gappy_site_data)
            self.init_config.update(init_params)

            # Update per-layer init for layered mode
            if self.spatial_mode == 'layered':
                self.layer_init = self.mapper.create_layered_initialization(gappy_site_data)

        # Initialize layer pools from site data (layered mode)
        if self.spatial_mode == 'layered':
            if gappy_site_data:
                self.layer_pools['ao_c'] = gappy_site_data.get('A0_c0', 5.0)
                self.layer_pools['ao_n'] = gappy_site_data.get('A0_n0', 0.2)
                self.layer_pools['sa_c'] = gappy_site_data.get('A_c0', 50.0)
                self.layer_pools['sa_n'] = gappy_site_data.get('A_n0', 2.5)
                self.layer_pools['sb_c'] = gappy_site_data.get('BL_c0', 20.0)
                self.layer_pools['sb_n'] = gappy_site_data.get('BL_n0', 1.0)
            else:
                self.layer_pools = {
                    'ao_c': 5.0, 'ao_n': 0.2,
                    'sa_c': 50.0, 'sa_n': 2.5,
                    'sb_c': 20.0, 'sb_n': 1.0,
                }

        # Create actual DEMENTpy Grid instances if enabled
        if self.enable_dement and DEMENT_AVAILABLE:
            print(f"Creating {self.n_grids} DEMENTpy Grid instance(s)...")

            # Grid configuration parameters
            gridsize = self.init_config.get('gridsize', 100)
            x = self.init_config.get('grid_x', 10)
            y = self.init_config.get('grid_y', 10)

            # Ensure gridsize matches x*y
            if gridsize != x * y:
                gridsize = x * y

            if self.spatial_mode == 'layered' and self.layer_init:
                # Create 3 grids with layer-specific conditions
                layer_names = ['AO', 'SA', 'SB']
                for i, layer_name in enumerate(layer_names):
                    try:
                        layer_cfg = self.layer_init[layer_name]
                        grid = DEMENTLibrary.create_grid(
                            end_time=365,
                            gridsize=gridsize, x=x, y=y,
                            n_taxa=self.init_config.get('n_taxa', 5),
                            n_substrates=self.init_config.get('n_substrates', 3),
                            substrate_c=layer_cfg['substrate_C'],
                            substrate_n=layer_cfg['substrate_N'],
                            substrate_p=self.init_config.get('litter_P', 5.0),
                            microbial_c=layer_cfg['microbial_C'],
                            temperature=layer_cfg.get('temperature', 15.0),
                            moisture_psi=self.init_config.get('moisture_psi', -0.5),
                            fb=self.init_config.get('fb', 0.1),
                            use_file_based_init=True
                        )
                        self.dement_grids.append(grid)
                        print(f"  Grid {i+1}/3 ({layer_name}) created")
                    except Exception as e:
                        print(f"  Error creating Grid {i+1} ({layer_name}): {e}")
                        print(f"    Falling back to empirical model")
                        self.enable_dement = False
                        break
            else:
                for i in range(self.n_grids):
                    try:
                        grid = DEMENTLibrary.create_grid(
                            end_time=365,
                            gridsize=gridsize, x=x, y=y,
                            n_taxa=self.init_config.get('n_taxa', 5),
                            n_substrates=self.init_config.get('n_substrates', 3),
                            substrate_c=self.init_config.get('litter_C', 500.0),
                            substrate_n=self.init_config.get('litter_N', 20.0),
                            substrate_p=self.init_config.get('litter_P', 5.0),
                            microbial_c=self.init_config.get('microbial_C', 50.0),
                            temperature=self.init_config.get('temperature', 15.0),
                            moisture_psi=self.init_config.get('moisture_psi', -0.5),
                            fb=self.init_config.get('fb', 0.1),
                            use_file_based_init=True
                        )
                        self.dement_grids.append(grid)
                        print(f"  Grid {i+1}/{self.n_grids} created")
                    except Exception as e:
                        print(f"  Error creating Grid {i+1}: {e}")
                        print(f"    Falling back to empirical model")
                        self.enable_dement = False
                        break
        else:
            # Create placeholder objects for fallback mode
            if self.spatial_mode == 'layered':
                layer_names = ['AO', 'SA', 'SB']
                layer_init = self.layer_init or self.mapper.DEFAULT_LAYER_INIT
                for i, layer_name in enumerate(layer_names):
                    layer_cfg = layer_init[layer_name]
                    grid_placeholder = {
                        'id': i,
                        'layer': layer_name,
                        'state': 'initialized',
                        'substrate_c': layer_cfg['substrate_C'],
                        'substrate_n': layer_cfg['substrate_N'],
                        'microbial_c': layer_cfg['microbial_C'],
                        'microbial_n': layer_cfg['microbial_N'],
                    }
                    self.dement_grids.append(grid_placeholder)
            else:
                for i in range(self.n_grids):
                    grid_placeholder = {
                        'id': i,
                        'state': 'initialized',
                        'substrate_c': self.init_config['litter_C'],
                        'substrate_n': self.init_config['litter_N'],
                        'microbial_c': self.init_config['microbial_C'],
                        'microbial_n': self.init_config['microbial_N'],
                    }
                    self.dement_grids.append(grid_placeholder)

        self.state_initialized = True
        if self.enable_dement:
            print(f"Initialized {self.n_grids} DEMENTpy Grid(s) - mechanistic mode")
        else:
            print(f"Initialized {self.n_grids} DEMENTpy Grid(s) - fallback mode")

    def couple_soil_decomp(self,
                          litter_c1: float, litter_c2: float,
                          litter_n1: float, litter_n2: float,
                          tempC: float, precip: float,
                          aow0_scaled_by_max: float,
                          saw0_scaled_by_fc: float,
                          sbw0_scaled_by_max: float,
                          plot_index: Optional[int] = None) -> Tuple[float, float]:
        """
        Main coupling function that replaces GAPPY's soil_decomp().

        This function matches the exact interface of soil_decomp() in soil.py,
        allowing it to be used as a drop-in replacement.

        Args:
            litter_c1: Aboveground litter C input (tc/ha/day)
            litter_c2: Belowground litter C input (tc/ha/day)
            litter_n1: Aboveground litter N input (tn/ha/day)
            litter_n2: Belowground litter N input (tn/ha/day)
            tempC: Daily temperature (°C)
            precip: Daily precipitation (cm/day)
            aow0_scaled_by_max: Scaled available water in AO layer
            saw0_scaled_by_fc: Scaled available water in SA layer
            sbw0_scaled_by_max: Scaled available water in SB layer
            plot_index: Plot index (0-199) for one_to_one mode

        Returns:
            Tuple of (avail_N, C_resp) where:
                avail_N: Available nitrogen (tn/ha)
                C_resp: CO2 respiration (tc/ha/day)
        """
        self.coupling_stats['total_calls'] += 1

        # Validate spatial mode requirements FIRST (before any processing)
        if self.spatial_mode == 'one_to_one':
            if plot_index is None:
                raise ValueError("plot_index required for one_to_one mode")
            if plot_index < 0 or plot_index >= self.n_plots:
                raise ValueError(f"plot_index must be between 0 and {self.n_plots-1}")

        # Initialize grids if not done yet
        if not self.state_initialized:
            self.initialize_dement_grids()

        # Increment day counter BEFORE processing (so fallback also tracks time)
        # This ensures proper time tracking regardless of which model is used
        current_day_for_run = self.current_day
        self.current_day += 1
        if self.current_day >= 365:
            self.current_day = 0
            self.current_year += 1

        # Route to layered decomposition if in layered mode
        if self.spatial_mode == 'layered':
            avail_N, C_resp = self._couple_layered_decomp(
                litter_c1, litter_c2, litter_n1, litter_n2,
                tempC, aow0_scaled_by_max, saw0_scaled_by_fc,
                sbw0_scaled_by_max
            )
            # Update statistics
            self.coupling_stats['total_litter_c'] += litter_c1 + litter_c2
            self.coupling_stats['total_litter_n'] += litter_n1 + litter_n2
            self.coupling_stats['total_resp'] += C_resp
            self.coupling_stats['total_n_avail'] += avail_N
            return avail_N, C_resp

        # If DEMENTpy not enabled, use simple fallback
        if not self.enable_dement:
            return self._fallback_soil_decomp(
                litter_c1, litter_c2, litter_n1, litter_n2,
                tempC, aow0_scaled_by_max, saw0_scaled_by_fc
            )

        # Convert inputs to DEMENTpy units
        total_litter_c = litter_c1 + litter_c2
        total_litter_n = litter_n1 + litter_n2

        # Update statistics
        self.coupling_stats['total_litter_c'] += total_litter_c
        self.coupling_stats['total_litter_n'] += total_litter_n

        # Select which grid(s) to run
        if self.spatial_mode == 'one_to_one':
            grid_indices = [plot_index]
        else:  # aggregated
            grid_indices = [0]

        # Run DEMENTpy for selected grids
        outputs = []
        for grid_idx in grid_indices:
            output = self._run_dement_daily_timestep(
                grid_idx, total_litter_c, total_litter_n,
                tempC, precip, aow0_scaled_by_max, saw0_scaled_by_fc,
                current_day_for_run
            )
            outputs.append(output)

        # Aggregate outputs if multiple grids
        if len(outputs) > 1:
            aggregated = self.converter.aggregate_plot_outputs(outputs)
        else:
            aggregated = outputs[0]

        # Update statistics
        self.coupling_stats['total_resp'] += aggregated['C_resp']
        self.coupling_stats['total_n_avail'] += aggregated['avail_N']

        return aggregated['avail_N'], aggregated['C_resp']

    def _run_dement_daily_timestep(self,
                                   grid_idx: int,
                                   daily_litter_c: float,
                                   daily_litter_n: float,
                                   daily_temp: float,
                                   daily_precip: float,
                                   daily_moisture_ao: float,
                                   daily_moisture_sa: float,
                                   day_of_year: int) -> Dict[str, float]:
        """
        Run DEMENTpy for one daily timestep.

        Runs actual DEMENTpy Grid timestep if enabled, otherwise uses fallback.

        Args:
            grid_idx: Index of the grid to run
            daily_litter_c: Daily litter C input (tc/ha/day)
            daily_litter_n: Daily litter N input (tn/ha/day)
            daily_temp: Daily temperature (°C)
            daily_precip: Daily precipitation (cm/day)
            daily_moisture_ao: Daily moisture in AO layer (scaled 0-1)
            daily_moisture_sa: Daily moisture in SA layer (scaled 0-1)
            day_of_year: Current day of year (0-364)

        Returns:
            Dictionary with 'avail_N' (tn/ha) and 'C_resp' (tc/ha/day)
        """
        # Convert to DEMENTpy units (g/m²/day)
        substrate = self.converter.litter_to_substrate(daily_litter_c, daily_litter_n)

        if self.enable_dement and DEMENT_AVAILABLE:
            # Run actual DEMENTpy mechanistic model
            try:
                grid = self.dement_grids[grid_idx]

                # Add daily litter inputs to grid substrates
                # Note: DEMENTpy Grid.Substrates is a DataFrame with C, N, P columns
                # We need to add the daily litter to the appropriate substrate pools
                # For now, add to 'DeadMic' pool which represents dead organic matter
                is_dead_mic = grid.Substrates.index == 'DeadMic'
                if is_dead_mic.any():
                    # Add litter C and N to the grid (distribute across grid cells)
                    n_cells = len(grid.Substrates.columns)
                    litter_c_per_cell = substrate['C'] / n_cells
                    litter_n_per_cell = substrate['N'] / n_cells

                    grid.Substrates.loc[is_dead_mic, 'C'] += litter_c_per_cell
                    grid.Substrates.loc[is_dead_mic, 'N'] += litter_n_per_cell

                # Run single daily timestep
                outputs = DEMENTLibrary.run_timestep(grid, day_of_year)

                # Extract daily outputs (already in g/m²)
                daily_n_avail_g_m2 = outputs['available_N']
                daily_resp_g_m2 = outputs['respiration']

                # Convert back to GAPPY units
                avail_N_tn_ha = self.converter.dement_to_gappy_n(daily_n_avail_g_m2)
                C_resp_tc_ha_day = self.converter.dement_to_gappy_resp(daily_resp_g_m2)

                return {
                    'avail_N': avail_N_tn_ha,
                    'C_resp': C_resp_tc_ha_day
                }

            except Exception as e:
                print(f"Warning: DEMENTpy execution failed: {e}")
                print("  Falling back to empirical model for this timestep")
                # Fall through to empirical model

        # Fallback: Simple empirical model similar to original GAPPY
        daily_substrate_c = substrate['C']  # g/m²/day
        daily_substrate_n = substrate['N']  # g/m²/day

        # Temperature adjustment
        if daily_temp >= -5.0:
            temp_factor = 3.0 ** (0.1 * (daily_temp - 1.0))
        else:
            temp_factor = 0.0

        # Moisture adjustment
        moisture_factor = max(0.2, 1.0 - (1.0 - daily_moisture_sa / 0.8)**2)

        # Respiration (rough estimate: 40% of inputs respired daily)
        resp_fraction = 0.4 * temp_factor * moisture_factor
        daily_resp_g_m2 = daily_substrate_c * resp_fraction

        # Available N (rough estimate: 10% of N inputs become available daily)
        n_avail_fraction = 0.1 * temp_factor * moisture_factor
        daily_n_avail_g_m2 = daily_substrate_n * n_avail_fraction

        # Convert back to GAPPY units
        avail_N_tn_ha = self.converter.dement_to_gappy_n(daily_n_avail_g_m2)
        C_resp_tc_ha_day = self.converter.dement_to_gappy_resp(daily_resp_g_m2)

        return {
            'avail_N': avail_N_tn_ha,
            'C_resp': C_resp_tc_ha_day
        }

    def _couple_layered_decomp(self,
                              litter_c1: float, litter_c2: float,
                              litter_n1: float, litter_n2: float,
                              tempC: float,
                              aow0_scaled_by_max: float,
                              saw0_scaled_by_fc: float,
                              sbw0_scaled_by_max: float) -> Tuple[float, float]:
        """
        Layered 3-grid decomposition: AO -> SA -> SB with cascading transfers.

        Each layer is processed by its own DEMENTpy grid (or fallback),
        with empirical transfer rates between layers matching GAPPY's
        original soil model structure.

        Returns:
            Tuple of (avail_N, C_resp) in GAPPY units (tn/ha, tc/ha/day)
        """
        ao_c = self.layer_pools['ao_c']
        ao_n = self.layer_pools['ao_n']
        sa_c = self.layer_pools['sa_c']
        sa_n = self.layer_pools['sa_n']
        sb_c = self.layer_pools['sb_c']
        sb_n = self.layer_pools['sb_n']

        # Temperature adjustment
        if tempC >= -5.0:
            tadjst = 3.0 ** (0.1 * (tempC - 1.0))
            tadjst1 = 2.5 ** (0.1 * (self._dampen_temperature(tempC, 0.8) - 1.0))
        else:
            tadjst = 0.0
            tadjst1 = 0.0

        # Moisture functions
        aow0_clamped = min(aow0_scaled_by_max, 0.5)
        aofunc = max((1.0 - (1.0 - aow0_clamped / 0.3) ** 2), 0.2)
        safunc = max(1.0 - (1.0 - saw0_scaled_by_fc / 0.8) ** 2, 0.2)

        # ---- AO layer: add aboveground litter, compute respiration ----
        ao_c += litter_c1
        ao_n += litter_n1
        ao_cn = ao_c / ao_n if ao_n > 0 else self.AO_CN_0

        if self.enable_dement and DEMENT_AVAILABLE and len(self.dement_grids) >= 3:
            # Run AO grid (grid 0)
            ao_output = self._run_dement_daily_timestep(
                0, litter_c1, litter_n1,
                tempC, 0.0, aow0_scaled_by_max, saw0_scaled_by_fc,
                self.current_day
            )
            resp1 = ao_output['C_resp']
        else:
            resp1 = tadjst * aofunc * self.AO_RESP * ao_c

        # Empirical AO -> SA transfer
        yxdn = resp1 / ao_cn if ao_cn > 0 else 0.0
        yxdc = yxdn * self.AO_CN_0

        # Cap total removal to avoid draining pool below floor
        total_ao_removal = yxdc + resp1
        if total_ao_removal > ao_c * 0.95:
            scale = (ao_c * 0.95) / total_ao_removal if total_ao_removal > 0 else 0.0
            resp1 *= scale
            yxdc *= scale
            yxdn *= scale

        ao_c = ao_c - yxdc - resp1
        ao_n = ao_n - yxdn

        # ---- SA layer: add belowground litter + AO transfer ----
        sa_c += yxdc + litter_c2
        sa_n += yxdn + litter_n2
        sa_cn = sa_c / sa_n if sa_n > 0 else self.SB_CN_0

        if self.enable_dement and DEMENT_AVAILABLE and len(self.dement_grids) >= 3:
            # Run SA grid (grid 1)
            sa_input_c = yxdc + litter_c2
            sa_input_n = yxdn + litter_n2
            sa_output = self._run_dement_daily_timestep(
                1, sa_input_c, sa_input_n,
                self._dampen_temperature(tempC, 0.8), 0.0,
                aow0_scaled_by_max, saw0_scaled_by_fc,
                self.current_day
            )
            resp2 = sa_output['C_resp']
        else:
            resp2 = tadjst1 * safunc * self.SA_RESP * sa_c

        # Empirical available N from SA
        avail_N = resp2 / sa_cn * max(0.5, (sa_cn - self.SA_CN_0) / sa_cn) if sa_cn > 0 else 0.0

        # Empirical SA -> SB transfer
        tosb = resp2 / self.SB_CN_0

        sa_c = sa_c - resp2 - tosb
        sa_n = sa_n - avail_N

        # ---- SB layer: add SA transfer, compute deep respiration ----
        sb_c += tosb

        if self.enable_dement and DEMENT_AVAILABLE and len(self.dement_grids) >= 3:
            # Run SB grid (grid 2)
            sb_output = self._run_dement_daily_timestep(
                2, tosb, tosb / self.SB_CN_0 if self.SB_CN_0 > 0 else 0.0,
                self._dampen_temperature(tempC, 0.6), 0.0,
                aow0_scaled_by_max, saw0_scaled_by_fc,
                self.current_day
            )
            resp3 = sb_output['C_resp']
        else:
            resp3 = sb_c * self.SB_RESP * tadjst1

        sb_c = sb_c - resp3

        # Total respiration
        C_resp = resp1 + resp2 + resp3

        # Clamp pools to safe minimums
        # AO needs a floor > 0 because GAPPY's soil_water uses A0_c0 to compute
        # water holding capacity (aow_max = A0_c0 * AO_MAX); zero causes div-by-zero
        ao_c = max(ao_c, 0.01)
        ao_n = max(ao_n, 1e-4)
        sa_c = max(sa_c, 0.01)
        sa_n = max(sa_n, 1e-4)
        sb_c = max(sb_c, 0.0)
        sb_n = max(sb_n, 1e-10)

        # Update layer pools
        self.layer_pools['ao_c'] = ao_c
        self.layer_pools['ao_n'] = ao_n
        self.layer_pools['sa_c'] = sa_c
        self.layer_pools['sa_n'] = sa_n
        self.layer_pools['sb_c'] = sb_c
        self.layer_pools['sb_n'] = sb_n

        return avail_N, C_resp

    @staticmethod
    def _dampen_temperature(surface_temp: float, depth_factor: float) -> float:
        """
        Attenuate temperature toward annual mean with depth.

        Deeper soil layers experience less temperature variability,
        trending toward an assumed annual mean of ~8 C.

        Args:
            surface_temp: Surface temperature (C)
            depth_factor: Attenuation factor (AO=1.0, SA=0.8, SB=0.6)

        Returns:
            Dampened temperature (C)
        """
        annual_mean = 8.0  # Assumed annual mean soil temperature
        return annual_mean + (surface_temp - annual_mean) * depth_factor

    def _fallback_soil_decomp(self,
                             litter_c1: float, litter_c2: float,
                             litter_n1: float, litter_n2: float,
                             tempC: float,
                             aow0_scaled_by_max: float,
                             saw0_scaled_by_fc: float) -> Tuple[float, float]:
        """
        Fallback to simple empirical soil decomposition model.

        This mimics the original GAPPY soil_decomp() for comparison.
        Used when enable_dement=False.

        Returns:
            Tuple of (avail_N, C_resp)
        """
        # Simple empirical constants (from original GAPPY)
        AO_RESP = 5.24e-4
        SA_RESP = 1.24e-5
        AO_CN_0 = 30.0
        SA_CN_0 = 4.0

        # Total inputs
        total_litter_c = litter_c1 + litter_c2
        total_litter_n = litter_n1 + litter_n2

        # Temperature adjustment
        if tempC >= -5.0:
            tadjst = 3.0**(0.1 * (tempC - 1.0))
        else:
            tadjst = 0.0

        # Moisture adjustment
        aofunc = max((1.0 - (1.0 - aow0_scaled_by_max / 0.3)**2), 0.2)
        safunc = max(1.0 - (1.0 - saw0_scaled_by_fc / 0.8)**2, 0.2)

        # Estimate soil C pools (simplified)
        ao_c = total_litter_c * 100.0  # Rough estimate
        sa_c = ao_c * 10.0

        # Respiration
        resp_ao = tadjst * aofunc * AO_RESP * ao_c
        resp_sa = tadjst * safunc * SA_RESP * sa_c
        total_resp = resp_ao + resp_sa

        # Available N (simplified)
        cn_ratio = total_litter_c / total_litter_n if total_litter_n > 0 else 30.0
        avail_n = resp_sa / cn_ratio * max(0.5, (cn_ratio - SA_CN_0) / cn_ratio)

        return avail_n, total_resp

    def get_state(self) -> Dict[str, Any]:
        """
        Get current state of all DEMENTpy grids.

        Returns:
            State dictionary that can be saved/restored
        """
        return {
            'year': self.current_year,
            'day': self.current_day,
            'grids': [grid for grid in self.dement_grids],
            'stats': self.coupling_stats.copy()
        }

    def set_state(self, state: Dict[str, Any]):
        """
        Restore DEMENTpy grids from saved state.

        Args:
            state: State dictionary from get_state()
        """
        self.current_year = state['year']
        self.current_day = state.get('day', 0)  # Default to 0 for old state files
        self.dement_grids = state['grids']
        self.coupling_stats = state['stats'].copy()
        self.state_initialized = True

    def advance_year(self):
        """Advance to next year and potentially reinitialize grids (pulse)."""
        self.current_year += 1

        # If using pulse structure, reinitialize periodically
        if self.runtime_config.get('pulse', 1) > 1:
            if self.current_year % self.runtime_config['pulse'] == 0:
                print(f"Pulse reinitialization at year {self.current_year}")
                # TODO: Reinitialize grids while preserving some state

    def print_statistics(self):
        """Print coupling statistics."""
        print("\n" + "="*60)
        print("DEMENTpy Coupling Statistics")
        print("="*60)
        print(f"Simulation time: Year {self.current_year}, Day {self.current_day}")
        print(f"Total coupling calls: {self.coupling_stats['total_calls']}")
        print(f"Total litter C input: {self.coupling_stats['total_litter_c']:.2f} tc/ha")
        print(f"Total litter N input: {self.coupling_stats['total_litter_n']:.2f} tn/ha")
        print(f"Total respiration: {self.coupling_stats['total_resp']:.2f} tc/ha")
        print(f"Total N available: {self.coupling_stats['total_n_avail']:.2f} tn/ha")

        if self.coupling_stats['total_litter_c'] > 0:
            resp_fraction = self.coupling_stats['total_resp'] / self.coupling_stats['total_litter_c']
            print(f"Respiration fraction: {resp_fraction:.2%}")

        if self.coupling_stats['total_litter_n'] > 0:
            n_avail_fraction = self.coupling_stats['total_n_avail'] / self.coupling_stats['total_litter_n']
            print(f"N availability fraction: {n_avail_fraction:.2%}")

        print("="*60 + "\n")


if __name__ == "__main__":
    # Demonstration of adapter usage
    print("DEMENTpy Adapter Demonstration\n")

    # Create adapter in aggregated mode
    adapter = DEMENTpyAdapter(
        n_plots=200,
        spatial_mode='aggregated',
        enable_dement=False  # Use fallback for now
    )

    # Simulate one coupling call
    litter_c1 = 0.01  # tc/ha/day aboveground
    litter_c2 = 0.005  # tc/ha/day belowground
    litter_n1 = 0.0005  # tn/ha/day
    litter_n2 = 0.00025  # tn/ha/day
    temp = 15.0  # °C
    precip = 0.3  # cm/day
    moisture_ao = 0.5
    moisture_sa = 0.7
    moisture_sb = 0.8

    print("Calling couple_soil_decomp() with:")
    print(f"  Litter C: {litter_c1 + litter_c2:.4f} tc/ha/day")
    print(f"  Litter N: {litter_n1 + litter_n2:.6f} tn/ha/day")
    print(f"  Temperature: {temp}°C")
    print(f"  Precipitation: {precip} cm/day\n")

    avail_n, resp = adapter.couple_soil_decomp(
        litter_c1, litter_c2, litter_n1, litter_n2,
        temp, precip, moisture_ao, moisture_sa, moisture_sb
    )

    print(f"Results:")
    print(f"  Available N: {avail_n:.6f} tn/ha")
    print(f"  Respiration: {resp:.6f} tc/ha/day")

    # Print statistics
    adapter.print_statistics()
