# Integration Module

This module provides the coupling layer between GAPPY (vegetation model) and DEMENTpy (microbial decomposition model) for the GoGs (Gap of Gaps) integrated ecosystem model.

## Overview

The integration layer replaces GAPPY's simple empirical soil decomposition function with DEMENTpy's mechanistic, enzyme-explicit microbial dynamics model.

## Components

### 1. `dement_adapter.py`
Main adapter class that orchestrates the coupling.

**Key class:** `DEMENTpyAdapter`
- Manages DEMENTpy Grid instances (1 or 200 depending on spatial mode)
- Provides `couple_soil_decomp()` method matching GAPPY's interface
- Handles temporal coupling (365 daily DEMENTpy steps per GAPPY year)
- Aggregates outputs from multiple grids

### 2. `unit_conversions.py`
Utilities for converting between GAPPY and DEMENTpy units.

**Key class:** `UnitConverter`
- GAPPY uses: tc/ha, tn/ha, cm/day
- DEMENTpy uses: g/m², g/m²/day
- Conversion factors: 1 tc/ha = 100 g/m²

### 3. `parameter_mapping.py`
Maps GAPPY parameters to DEMENTpy initialization.

**Key class:** `ParameterMapper`
- Creates DEMENTpy runtime configurations
- Maps GAPPY soil C/N pools to DEMENTpy initial conditions
- Provides default forest soil parameters

## Usage

### Basic Usage (Aggregated Mode)

```python
from model.integration import DEMENTpyAdapter

# Create adapter in aggregated mode (all plots → 1 DEMENTpy grid)
adapter = DEMENTpyAdapter(
    n_plots=200,
    spatial_mode='aggregated',
    enable_dement=False  # Use fallback empirical model for now
)

# Call from GAPPY's soil decomposition loop
avail_n, resp = adapter.couple_soil_decomp(
    litter_c1=0.01,    # Aboveground litter C (tc/ha/day)
    litter_c2=0.005,   # Belowground litter C (tc/ha/day)
    litter_n1=0.001,   # Aboveground litter N (tn/ha/day)
    litter_n2=0.0005,  # Belowground litter N (tn/ha/day)
    tempC=15.0,        # Temperature (°C)
    precip=0.3,        # Precipitation (cm/day)
    aow0_scaled_by_max=0.5,  # AO layer moisture
    saw0_scaled_by_fc=0.7,   # SA layer moisture
    sbw0_scaled_by_max=0.8   # SB layer moisture
)
```

### One-to-One Mode (Spatially Explicit)

```python
# Create adapter with one DEMENTpy grid per GAPPY plot
adapter = DEMENTpyAdapter(
    n_plots=200,
    spatial_mode='one_to_one',
    enable_dement=False
)

# Must provide plot_index in one-to-one mode
for plot_idx in range(200):
    avail_n, resp = adapter.couple_soil_decomp(
        litter_c1, litter_c2, litter_n1, litter_n2,
        tempC, precip, moisture_ao, moisture_sa, moisture_sb,
        plot_index=plot_idx  # Specify which plot/grid
    )
```

### Custom Configuration

```python
# Custom DEMENTpy runtime parameters
runtime_config = {
    'pulse': 1,
    'end_time': 365,
    'interval': 30,
    'grid_size': 10,
    'n_taxa': 5
}

# Custom initialization parameters
init_config = {
    'litter_C': 500.0,  # g/m²
    'litter_N': 20.0,
    'SOM_C': 5000.0,
    'temperature': 15.0,
    'moisture': 0.6
}

adapter = DEMENTpyAdapter(
    spatial_mode='aggregated',
    runtime_config=runtime_config,
    init_config=init_config,
    enable_dement=True  # Use actual DEMENTpy when ready
)
```

## Integration with GAPPY

### Modifying soil.py

To integrate with GAPPY, modify `model/vegetation/soil.py`:

```python
from model.integration import DEMENTpyAdapter

class SoilData:
    def __init__(self, use_dement=False):
        # ... existing init ...

        self.use_dement = use_dement
        if use_dement:
            self.dement_adapter = DEMENTpyAdapter(
                n_plots=200,
                spatial_mode='aggregated',
                enable_dement=True
            )

    def soil_decomp(self, litter_c1, litter_c2, ...):
        if self.use_dement:
            return self.dement_adapter.couple_soil_decomp(
                litter_c1, litter_c2, litter_n1, litter_n2,
                tempC, precip, aow0_scaled_by_max,
                saw0_scaled_by_fc, sbw0_scaled_by_max
            )
        else:
            # Original empirical model
            ...
```

## Unit Conversions Reference

| GAPPY Unit | DEMENTpy Unit | Conversion |
|------------|---------------|------------|
| tc/ha | g C/m² | × 100 |
| tn/ha | g N/m² | × 100 |
| cm/day | cm/day | direct |
| °C | °C | direct |

## Current Status

**Phase 1 (Complete):** Adapter layer implementation
- ✅ Directory structure created
- ✅ `DEMENTpyAdapter` class implemented
- ✅ Unit conversion utilities
- ✅ Parameter mapping system
- ✅ Fallback empirical model

**Phase 2 (Pending):** Modify soil_decomp() in GAPPY
**Phase 3 (Pending):** Modify DEMENTpy for library usage
**Phase 4 (Pending):** Testing and validation

## Testing

Run the adapter demonstration:

```bash
cd model/integration
python dement_adapter.py
```

Run unit conversion tests:

```bash
python unit_conversions.py
```

Run parameter mapping demonstration:

```bash
python parameter_mapping.py
```

## Known Limitations

1. **DEMENTpy not yet integrated:** Currently uses fallback empirical model
2. **No true spatial heterogeneity:** Even in one-to-one mode, grids are independent
3. **Annual timestep:** DEMENTpy runs 365 days per GAPPY year (could be optimized)
4. **Memory usage:** One-to-one mode with 200 grids may be memory-intensive

## Future Enhancements

1. Implement actual DEMENTpy Grid initialization and execution
2. Add parallel processing for one-to-one mode
3. Implement adaptive timestep coupling
4. Add comprehensive mass balance checking
5. Create visualization tools for coupled outputs
6. Optimize performance with Cython/numba

## References

- See `INTEGRATION_PLAN.md` for complete integration strategy
- GAPPY documentation: `model/vegetation/`
- DEMENTpy documentation: `model/microbiome/DEMENTpy/`

---

**Version:** 0.1.0 (Phase 1 Complete)
**Status:** Adapter layer ready, awaiting DEMENTpy library interface
