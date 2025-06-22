# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GoGs (Gap of Gaps) is a predictive synthesis framework that models vegetation and microbiome interactions in ecosystems and the biosphere. The project aims to create an integrative ecosystem model using a Pattern, Diversity, and Process (PDP) framework.

## Architecture

The codebase consists of three main components:

### 1. Vegetation Model (Fortran-based)
- Located in `model/vegetation/src/`
- Based on UVAFME (University of Virginia Forest Model Enhanced)
- Core Fortran modules: Climate, Species, Tree, Site, Soil, Model, etc.
- Entry point: `UVAFME.f90`
- Python interface built using f2py and Meson

### 2. Microbiome Model (Python-based)
- Will be based on DEMENTpy
- Currently in development

### 3. Integration Layer (Python)
- Main entry points: `main.py` and `model/main.py`
- Uses numpy for numerical operations
- Python bindings to Fortran vegetation model

## Building and Compilation

### Vegetation Model (Fortran)
The Fortran vegetation model can be built using multiple approaches:

**Option 1: Traditional Makefile**
```bash
cd model/vegetation/src/
make
```

**Option 2: Python f2py with Meson backend**
```bash
cd model/vegetation/
python compile.py
```

**Option 3: Meson build system**
```bash
cd model/vegetation/src/
meson setup builddir
meson compile -C builddir
```

### Python Environment
```bash
# Install dependencies using uv (recommended)
uv sync

# Or using pip
pip install -e .
```

## Running the Model

### Basic execution
```bash
python main.py
```

### Vegetation model integration
```bash
cd model/
python main.py
```

## Dependencies

- Python >=3.9
- numpy >1.25.2
- meson >=1.7.2
- Fortran compiler (ifort recommended, gfortran also supported)
- f2py (part of numpy)

## File Structure

- `model/vegetation/input_data/`: Climate and site data (CSV files)
- `model/vegetation/output_data/`: Model outputs (Species, Genus, Carbon/Nitrogen data)
- `model/vegetation/src/`: Fortran source code
- `manuscript/`: Documentation and drafts
- `reaction_diffusion_model.ipynb`: Jupyter notebook for model development

## Development Notes

- The vegetation model is a port/adaptation of UVAFME
- Fortran-Python interface uses f2py for wrapping
- No formal test suite currently exists
- Input data includes climate, species, and site information
- Model outputs ecosystem dynamics over time periods specified in runtime configuration