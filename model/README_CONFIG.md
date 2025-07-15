# UVAFME Configuration Guide

## Configuration File: `uvafme_config.json`

The UVAFME model uses a JSON configuration file to set simulation parameters. The model looks for this file in the following locations:

1. `uvafme_config.json` (current directory)
2. `input_data/uvafme_config.json` 

If no configuration file is found, the model will use default parameters.

## Configuration Parameters

### Basic Simulation Parameters
- `numyears`: Number of years to simulate (default: 50)
- `numplots`: Number of plots per site (default: 1)
- `maxtrees`: Maximum number of trees per plot (default: 1000)
- `maxheight`: Maximum tree height in meters (default: 50)
- `year_print_interval`: How often to write output data in years (default: 10)

### Random Number Generation
- `fixed_seed`: Use fixed random seed for reproducible results (default: true)
- `same_climate`: Use same climate for all runs (default: true)
- `debug`: Enable debug output (default: false)

### Spinup Configuration
- `spinup`: Enable spinup period (default: false)
- `spinup_yrs`: Number of spinup years (default: 20)

### Plot/Site Parameters
- `plotsize`: Plot size in hectares (default: 1.0)
- `rootdepth`: Root depth in meters (default: 2.0)
- `fire_level`: Fire probability (default: 0.01)
- `wind_level`: Wind probability (default: 0.05)

### Output Options
- `tree_level_data`: Write individual tree data (default: false) **Warning: Large files!**
- `plot_level_data`: Write plot-level data (default: false)

### Climate Change (Optional)
- `with_clim_change`: Enable climate change scenario (default: false)
- `linear_cc`: Use linear climate change (default: false)
- `use_gcm`: Use GCM data (default: false)
- `begin_change_year`: Year to start climate change (default: 30)
- `duration_of_change`: Duration of change in years (default: 20)
- `incr_or_decr`: "incr" for increase, "decr" for decrease (default: "incr")

### Climate Change Amounts
- `incr_tmin_by`: Increase minimum temperature by °C (default: 2.0)
- `incr_tmax_by`: Increase maximum temperature by °C (default: 2.0)
- `incr_precip_by`: Increase precipitation by fraction (default: 0.1)
- `decr_tmin_by`: Decrease minimum temperature by °C (default: 0.0)
- `decr_tmax_by`: Decrease maximum temperature by °C (default: 0.0)
- `decr_precip_by`: Decrease precipitation by fraction (default: 0.0)

## Example Configuration

```json
{
  "numyears": 100,
  "numplots": 2,
  "maxtrees": 500,
  "maxheight": 45,
  "year_print_interval": 5,
  "fixed_seed": true,
  "spinup": true,
  "spinup_yrs": 30,
  "plotsize": 0.5,
  "with_clim_change": true,
  "begin_change_year": 50,
  "incr_tmin_by": 3.0,
  "incr_tmax_by": 3.0,
  "incr_precip_by": 0.2
}
```

## Creating Configuration Files

### Method 1: Use the utility script
```bash
uv run python create_config.py
```

### Method 2: Python code
```python
from vegetation import params
params.create_default_config("my_config.json")
```

### Method 3: Load from dictionary
```python
from vegetation import params

config = {
    "numyears": 75,
    "maxtrees": 800,
    "year_print_interval": 5
}

params.load_from_dict(config)
```

## Running with Configuration

The model automatically loads the configuration file:

```bash
uv run python main.py
```

## Output Files

The model creates several output files in the `output_data/` directory:

- `species_data.csv`: Species-level biomass and tree counts
- `site_data.csv`: Site-level climate and water data
- `soil_data.csv`: Soil carbon and nitrogen data
- `trees_year_X.csv`: Individual tree data (if enabled)

## Tips

1. **Start small**: Use shorter simulations (20-50 years) for testing
2. **Tree data warning**: Enabling `tree_level_data` creates very large files
3. **Reproducibility**: Keep `fixed_seed: true` for consistent results
4. **Climate change**: Enable gradually - start with small changes
5. **Performance**: Reduce `maxtrees` if simulation is slow

## Validation

The model validates configuration parameters on startup. Common issues:
- Negative values for years, trees, or heights
- Print intervals larger than simulation years
- Invalid climate change parameters