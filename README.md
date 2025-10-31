# GoGs: Gap of Gaps
**A predictive synthesis framework of vegetation and microbiome underlying ecosystem and the Biosphere functioning and dynamics.**

Theoretical and empirical research have established a view of ecosystems and the biosphere as complex adaptive systems since the early 1990s. After half of a century of effort in modelling ecosystems and the biosphere interacting with the atmosphere, a premise remains of up-scaling biogeochemical processes without appreciating their complexity to capture self-organization and emergence. However, this prevailing and fundamental assumption, however feasible and impactful, is losing its grounds as anthropogenic disturbances facing ecosystems and the biosphere become more frequent and intense. We argue that the time comes of formulating a real integrative ecosystem model that can transcend boundaries of scale and organism (especially plant and microbe, two major branches of the tree of life), explain their high-dimensional interactions with the environment, and ideally march together with empirical research. Here, we propose a framework of Pattern, Diversity, and Process (PDP) for demography-based integrative modelling of ecosystems and the biosphere to contribute to this grand goal.

---

## Architecture

The GoGs framework integrates two major model components:

### 1. **GAPPY** (Vegetation Model)
- Based on [GAPpy](https://github.com/bioatmosphere/GAPpy.git)
- Located in `model/vegetation/`
- Simulates forest gap dynamics, tree growth, species competition, and succession
- Provides litter inputs (leaf, woody debris) to the microbiome model

### 2. **DEMENTpy** (Microbiome Model)
- Based on [DEMENTpy](https://github.com/DEMENT-Model/DEMENTpy)
- Located in `model/microbiome/DEMENTpy/`
- Simulates microbial communities, substrate degradation, and nutrient cycling
- Generates mineral nutrients (NH4, PO4) for vegetation uptake

### 3. **Integration Layer**
- Located in `model/integration/`
- Bidirectional coupling between vegetation and microbiome
- Unit conversions and parameter mapping
- Manages data exchange between models

## Project Structure

```
GoG/
├── model/
│   ├── vegetation/          # GAPPY vegetation model (Python)
│   │   ├── gappy.py        # Main model class
│   │   ├── species.py      # Species data and traits
│   │   ├── site.py         # Site conditions
│   │   ├── soil.py         # Soil processes
│   │   ├── tree.py         # Individual tree dynamics
│   │   └── src/            # Original Fortran source (UVAFME)
│   ├── microbiome/         # DEMENTpy microbiome model
│   │   └── DEMENTpy/       # Submodule
│   └── integration/        # Coupling layer
│       ├── parameter_mapping.py
│       └── unit_conversions.py
├── input_data/             # Model input files
├── output_data/            # Model outputs
├── run_with_config.py      # Run with configuration file
└── run_integrated_model.py # Run fully integrated model
```

## Installation

### Requirements
- Python ≥3.9
- numpy ≥1.25.2
- pandas
- matplotlib (for plotting)

### Setup

```bash
# Clone the repository
git clone <repository-url>
cd GoG

# Install dependencies using uv (recommended)
uv sync

# Or using pip
pip install -e .
```

## Usage

### Quick Start

```bash
# Run with default configuration
uv run python run_with_config.py

# Run integrated vegetation-microbiome model
uv run python run_integrated_model.py
```

### Configuration

Model parameters are configured through JSON files:

- **Vegetation**: `input_data/gappy_config.json`
- **Microbiome**: Set programmatically in run scripts

Example configuration:

```python
from model.vegetation import GAPPYModel

model = GAPPYModel()
model.initialize_input_files()
model.run_model()
```

## Key Features

### Vegetation Model (GAPPY)
- Individual-based forest gap model
- Species-specific traits (growth rates, shade tolerance, drought tolerance)
- Climate-driven dynamics (temperature, precipitation)
- Produces organic matter inputs (litter, woody debris)

### Microbiome Model (DEMENTpy)
- Spatially explicit microbial community dynamics
- Enzyme-mediated substrate degradation
- Organic and mineral monomer cycling
- Emergent carbon use efficiency (CUE)

### Integration
- **Vegetation → Microbiome**: Litter inputs (leaf and woody material with C:N:P ratios)
- **Microbiome → Vegetation**: Mineral nutrients (NH4, PO4) from decomposition
- Unit conversions: g C/m²/yr ↔ daily rates
- Temporal coupling: Annual vegetation → daily microbiome → annual feedback

## Monomer Types

The model distinguishes between:

- **Mineral Monomers** (NH4, PO4): Inorganic nutrients with fixed stoichiometry
  - Generated through microbial mineralization
  - Subject to leaching losses
  - Pure N or P (no carbon)

- **Organic Monomers**: Complex compounds from substrate degradation
  - Variable C:N:P ratios inherited from parent substrates
  - Contain carbon (energy source)
  - Produced by enzymatic depolymerization

## Output

Model outputs include:
- Species composition and biomass over time
- Soil carbon and nitrogen pools
- Microbial community dynamics
- Nutrient fluxes
- Ecosystem respiration and carbon use efficiency

## Development Status

This is an active research project integrating vegetation and microbiome models. The codebase includes:
- ✅ Python translation of UVAFME → GAPPY
- ✅ DEMENTpy microbiome model integration
- ✅ Bidirectional coupling layer
- 🚧 Full validation and calibration ongoing

## References

- **GAPPY**: [GitHub Repository](https://github.com/bioatmosphere/GAPpy.git)
- **DEMENTpy**: [GitHub Repository](https://github.com/DEMENT-Model/DEMENTpy)

## Citation

[Citation information to be added]

## License

[License information to be added]
