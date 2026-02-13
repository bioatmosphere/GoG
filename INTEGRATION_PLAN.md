# GAPPY-DEMENTpy Integration Plan
## GoGs (Gap of Gaps) Model Integration Strategy

**Version:** 1.0
**Date:** 2025-01-28
**Status:** Planning Phase

---

## Executive Summary

This document outlines the integration strategy for coupling the GAPPY forest gap model (vegetation dynamics) with DEMENTpy (microbial decomposition dynamics) to create the GoGs (Gap of Gaps) integrated ecosystem model.

### Key Goals
1. Replace GAPPY's simple empirical soil decomposition with mechanistic DEMENTpy model
2. Maintain bidirectional coupling: vegetation → microbiome (litter inputs) and microbiome → vegetation (nutrient availability)
3. Ensure temporal and spatial scale compatibility
4. Preserve modularity for independent model testing and validation

---

## 1. Model Architecture Overview

### Current State

#### GAPPY (Vegetation Model)
- **Language:** Python (translated from Fortran)
- **Structure:** Forest gap model with individual tree dynamics
- **Spatial scale:** 200 plots × 500 m² each (1 hectare total)
- **Temporal scale:** Yearly timesteps (365 daily cycles internally)
- **Key outputs for coupling:**
  - Litter C and N inputs (aboveground and belowground)
  - Temperature and precipitation
  - Soil water status
- **Key inputs from coupling:**
  - Available N for plant uptake
  - Soil respiration (CO₂)

#### DEMENTpy (Microbiome Model)
- **Language:** Python
- **Structure:** Grid-based microbial community model
- **Spatial scale:** Configurable grid (e.g., 10×10 cells)
- **Temporal scale:** Daily timesteps with pulse structure
- **Key outputs for coupling:**
  - Available N from decomposition
  - CO₂ respiration
  - Microbial biomass
- **Key inputs from coupling:**
  - Substrate (litter) C and N inputs
  - Temperature and moisture

---

## 2. Integration Strategy

### 2.1 Coupling Architecture

We propose a **hierarchical replacement coupling** approach:

```
┌─────────────────────────────────────────────────────────┐
│                    GAPPYModel                          │
│  ┌──────────────────────────────────────────────────┐  │
│  │         ForestModel (yearly timestep)            │  │
│  │  • Tree growth, competition, mortality           │  │
│  │  • Litter production                             │  │
│  │  • Nutrient uptake                               │  │
│  └─────────────────┬────────────────────────────────┘  │
│                    │                                    │
│  ┌─────────────────▼────────────────────────────────┐  │
│  │              SoilData                            │  │
│  │  ┌────────────────────────────────────────────┐ │  │
│  │  │  soil_decomp() [REPLACEMENT POINT]         │ │  │
│  │  │  Currently: Simple empirical model         │ │  │
│  │  │  Future: DEMENTpy coupling layer           │ │  │
│  │  └────────────────────────────────────────────┘ │  │
│  └──────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────┐
│                   DEMENTpyAdapter                       │
│  • Manages 200 DEMENTpy Grid instances                 │
│  • Spatial mapping: 1 forest plot → 1 microbial grid   │
│  • Temporal coupling: 1 year GAPPY → 365 days DEMENT  │
│  • Input/output transformation and aggregation         │
└─────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────┐
│               DEMENTpy (modified)                       │
│  • Grid-based microbial community dynamics             │
│  • Substrate degradation and enzyme production         │
│  • Nutrient cycling and CO₂ respiration                │
└─────────────────────────────────────────────────────────┘
```

### 2.2 Spatial Coupling Strategy

**Option A: One-to-One Mapping (Recommended)**
- Each GAPPY plot (500 m²) → 1 DEMENTpy grid instance
- 200 plots → 200 independent DEMENTpy simulations
- Advantages: Spatially explicit, captures plot heterogeneity
- Disadvantages: Higher computational cost

**Option B: Aggregated Approach**
- All 200 plots → 1 averaged DEMENTpy simulation
- Aggregate litter inputs across plots
- Advantages: Computationally efficient
- Disadvantages: Loss of spatial heterogeneity

**Recommendation:** Start with Option B for testing, transition to Option A for production

### 2.3 Temporal Coupling Strategy

**Timestep synchronization:**

```
GAPPY Year Loop (500 years)
├── Year N (365 daily cycles)
│   ├── Accumulate daily litter production
│   ├── Calculate daily climate averages
│   └── Day 365: Call soil_decomp()
│       └── DEMENTpy coupling point
│           ├── Convert accumulated annual litter → daily inputs
│           ├── Run DEMENTpy for 365 days
│           ├── Aggregate DEMENTpy outputs
│           └── Return: available_N, CO2_resp
│
└── Update soil C/N pools, continue to Year N+1
```

**Key decisions:**
1. **Frequency:** Call DEMENTpy once per GAPPY year (not daily) to balance accuracy and performance
2. **Litter distribution:** Distribute annual litter evenly across 365 days for DEMENTpy
3. **Climate forcing:** Use GAPPY's daily climate data to drive DEMENTpy temperature/moisture

---

## 3. Implementation Roadmap

### Phase 1: Adapter Layer Development (2-3 weeks)

**Deliverable:** `DEMENTpyAdapter` class

**Tasks:**
1. Create `model/integration/dement_adapter.py`
2. Implement DEMENTpy wrapper class:
   - Initialize DEMENTpy Grid with GAPPY-derived parameters
   - Map GAPPY litter inputs → DEMENTpy substrate inputs
   - Map DEMENTpy outputs → GAPPY soil variables
3. Handle unit conversions:
   - Litter: tc/ha/year → substrate pool units
   - N availability: DEMENTpy units → tn/ha
   - Respiration: DEMENTpy units → tc/ha/day

**Code structure:**
```python
class DEMENTpyAdapter:
    def __init__(self, runtime_config, spatial_mode='aggregated'):
        """
        Initialize adapter with configuration.

        Args:
            runtime_config: Dict with DEMENTpy runtime parameters
            spatial_mode: 'aggregated' or 'one_to_one'
        """
        self.spatial_mode = spatial_mode
        self.dement_grids = []  # List of Grid instances

        if spatial_mode == 'one_to_one':
            # Create 200 Grid instances (one per plot)
            for i in range(200):
                self.dement_grids.append(Grid(runtime_config, data_init))
        else:
            # Create single aggregated Grid
            self.dement_grids.append(Grid(runtime_config, data_init))

    def couple_soil_decomp(self, litter_c1, litter_c2, litter_n1, litter_n2,
                          tempC, precip, moisture_state, plot_index=None):
        """
        Replace soil_decomp() with DEMENTpy coupling.

        Matches GAPPY soil_decomp() interface exactly.

        Returns:
            avail_N: Available N (tn/ha)
            C_resp: CO2 respiration (tc/ha/day)
        """
        # Convert inputs
        substrates = self._litter_to_substrates(litter_c1, litter_c2,
                                                litter_n1, litter_n2)

        # Select appropriate Grid instance
        if self.spatial_mode == 'one_to_one':
            grid = self.dement_grids[plot_index]
        else:
            grid = self.dement_grids[0]

        # Run DEMENTpy for 365 days
        for day in range(365):
            grid.degradation(day)
            grid.uptake(day)
            grid.metabolism(day)
            grid.mortality(day)
            grid.reproduction(day)

        # Extract and convert outputs
        avail_N = self._extract_available_N(grid)
        C_resp = self._extract_respiration(grid)

        return avail_N, C_resp
```

### Phase 2: Modified soil_decomp() (1 week)

**Deliverable:** Updated `model/vegetation/soil.py`

**Tasks:**
1. Add adapter initialization to `SoilData.__init__()`
2. Create toggle between empirical and mechanistic models:
   ```python
   def __init__(self, use_dement=False):
       self.use_dement = use_dement
       if use_dement:
           self.dement_adapter = DEMENTpyAdapter(...)
   ```
3. Modify `soil_decomp()` to dispatch to adapter:
   ```python
   def soil_decomp(self, litter_c1, litter_c2, ...):
       if self.use_dement:
           return self.dement_adapter.couple_soil_decomp(
               litter_c1, litter_c2, litter_n1, litter_n2,
               tempC, precip, moisture_state, plot_index
           )
       else:
           # Original empirical model
           ...
   ```

### Phase 3: DEMENTpy Modifications (1-2 weeks)

**Deliverable:** Modified `model/microbiome/DEMENTpy/src/`

**Required changes:**
1. **Remove command-line dependency:**
   - Make Grid initialization accept parameters directly (not from sys.argv)
   - Remove file I/O requirements for initialization

2. **Create library interface:**
   ```python
   # New file: model/microbiome/DEMENTpy/src/coupling.py
   class DEMENTLibrary:
       """Library interface for coupled simulations."""

       @staticmethod
       def create_grid(params_dict):
           """Create Grid from parameter dictionary."""
           return Grid(runtime=params_dict, data_init=...)

       @staticmethod
       def get_outputs(grid):
           """Extract coupling-relevant outputs."""
           return {
               'available_N': ...,
               'CO2_respiration': ...,
               'microbial_biomass': ...
           }
   ```

3. **Add state persistence:**
   - Allow Grid to persist state between calls
   - Support annual reinitialization (pulse structure)

### Phase 4: Configuration and Parameters (1 week)

**Deliverable:** Parameter mapping framework

**Tasks:**
1. Create `model/integration/parameter_mapping.py`:
   - Map GAPPY site parameters → DEMENTpy initialization
   - Create default DEMENTpy runtime config for forest ecosystems

2. Add configuration options to `input_data/gappy_config.json`:
   ```json
   {
     "use_mechanistic_soil": true,
     "dement_config": {
       "spatial_mode": "aggregated",
       "pulse": 1,
       "end_time": 365,
       "interval": 30,
       "dispersal": 0
     }
   }
   ```

### Phase 5: Testing and Validation (2-3 weeks)

**Deliverable:** Validated integrated model

**Test scenarios:**
1. **Unit tests:**
   - Adapter input/output transformations
   - Unit conversions
   - Spatial mapping correctness

2. **Comparison tests:**
   - Run GAPPY with empirical soil model (baseline)
   - Run GAPPY with DEMENTpy (integrated)
   - Compare: Available N, respiration, biomass dynamics

3. **Sensitivity analysis:**
   - Test with different litter inputs
   - Test with different climate scenarios
   - Verify mass balance (C and N conservation)

4. **Performance testing:**
   - Measure computational overhead
   - Optimize bottlenecks
   - Profile both spatial modes

---

## 4. Technical Specifications

### 4.1 Interface Contract

**Input transformation (GAPPY → DEMENTpy):**

| GAPPY Variable | Units | DEMENTpy Variable | Units | Conversion |
|----------------|-------|-------------------|-------|------------|
| litter_c1 | tc/ha/day | Substrate C (above) | g C/m² | × 0.1 |
| litter_c2 | tc/ha/day | Substrate C (below) | g C/m² | × 0.1 |
| litter_n1 | tn/ha/day | Substrate N (above) | g N/m² | × 0.1 |
| litter_n2 | tn/ha/day | Substrate N (below) | g N/m² | × 0.1 |
| tempC | °C | Temperature | °C | direct |
| precip | cm/day | Moisture input | cm/day | direct |

**Output transformation (DEMENTpy → GAPPY):**

| DEMENTpy Variable | Units | GAPPY Variable | Units | Conversion |
|------------------|-------|----------------|-------|------------|
| Available N | g N/m² | avail_N | tn/ha | × 10 |
| CO₂ flux | g C/m²/day | C_resp | tc/ha/day | × 10 |

### 4.2 State Management

**DEMENTpy state persistence:**
- Microbial biomass carries over year-to-year
- Substrate pools update continuously
- Annual pulse reinitialization optional (configurable)

**GAPPY-DEMENTpy synchronization:**
```python
# Year N
gappy_state = {
    'A0_c0': ...,  # Soil C pools
    'A_c0': ...,
    'microbial_state': dement_adapter.get_state()  # NEW
}

# Year N+1
dement_adapter.set_state(gappy_state['microbial_state'])
```

### 4.3 Mass Balance Verification

**Conservation checks:**
1. Total ecosystem C = vegetation_C + soil_C + microbial_C + respiration
2. Total ecosystem N = vegetation_N + soil_N + microbial_N + losses
3. Implement mass balance logging in adapter

---

## 5. Expected Outcomes

### Scientific Benefits
1. **Mechanistic soil processes:** Replace empirical decomposition with enzyme-explicit microbial dynamics
2. **Microbial diversity effects:** Explore how microbial community composition affects forest dynamics
3. **Climate change responses:** Better predict ecosystem C and N cycling under warming
4. **Nutrient limitation:** More realistic representation of N availability for plants

### Model Improvements
1. **Bidirectional feedbacks:** Capture vegetation → microbiome → vegetation cycles
2. **Spatial heterogeneity:** Track plot-level microbial communities
3. **Temporal dynamics:** Daily microbial processes within annual forest dynamics
4. **Emergent properties:** Discover new ecosystem behaviors from coupled processes

---

## 6. Risks and Mitigation

### Technical Risks

**Risk 1: Computational Performance**
- **Impact:** 200 DEMENTpy instances may be too slow
- **Mitigation:**
  - Start with aggregated mode
  - Implement parallel processing for one-to-one mode
  - Profile and optimize bottlenecks
  - Consider GPU acceleration for DEMENTpy

**Risk 2: Numerical Instability**
- **Impact:** Coupling different timesteps may cause instability
- **Mitigation:**
  - Implement mass balance checks
  - Add numerical damping for litter inputs
  - Test with small timesteps first
  - Validate against known benchmarks

**Risk 3: Parameter Uncertainty**
- **Impact:** DEMENTpy initialization may not match forest soils
- **Mitigation:**
  - Literature review for forest soil parameters
  - Calibration against GAPPY default behavior
  - Sensitivity analysis for key parameters
  - Document parameter provenance

### Scientific Risks

**Risk 4: Model Mismatch**
- **Impact:** Spatial/temporal scales may be incompatible
- **Mitigation:**
  - Flexible coupling frequency (annual, seasonal, monthly)
  - Test multiple aggregation schemes
  - Compare with decoupled runs

**Risk 5: Validation Challenges**
- **Impact:** Difficult to validate coupled model without field data
- **Mitigation:**
  - Compare trends with long-term forest plots
  - Use literature C:N ratios as benchmarks
  - Test against GAPPY historical validation

---

## 7. Timeline

| Phase | Duration | Dependencies | Deliverables |
|-------|----------|--------------|--------------|
| 1. Adapter Layer | 2-3 weeks | DEMENTpy understanding | `DEMENTpyAdapter` class |
| 2. Modify soil_decomp | 1 week | Phase 1 | Updated `soil.py` |
| 3. DEMENTpy mods | 1-2 weeks | Phase 1 | Library interface |
| 4. Configuration | 1 week | Phase 2, 3 | Parameter files |
| 5. Testing | 2-3 weeks | Phase 4 | Validated model |
| **Total** | **7-10 weeks** | | Integrated GoGs model |

---

## 8. Success Criteria

### Milestone 1: Successful Compilation
- [ ] Adapter compiles without errors
- [ ] GAPPY runs with adapter (even if results incorrect)
- [ ] No runtime crashes

### Milestone 2: Mass Balance
- [ ] C and N conservation verified
- [ ] No negative pools
- [ ] Respiration within reasonable bounds

### Milestone 3: Validation
- [ ] Integrated model produces stable 500-year runs
- [ ] Available N dynamics are realistic
- [ ] Forest biomass trends comparable to empirical model

### Milestone 4: Scientific Output
- [ ] Model demonstrates emergent microbiome-vegetation feedbacks
- [ ] Produces novel predictions not possible with empirical model
- [ ] Ready for scientific publication

---

## 9. Next Steps

### Immediate Actions (Week 1)
1. Create `model/integration/` directory structure
2. Draft `DEMENTpyAdapter` skeleton class
3. Identify DEMENTpy initialization requirements
4. Set up test framework for adapter

### Short-term Goals (Month 1)
1. Complete Phase 1 (Adapter Layer)
2. Implement basic coupling with aggregated mode
3. Run first integrated test simulation
4. Document initial results and issues

### Long-term Goals (Months 2-3)
1. Complete all implementation phases
2. Comprehensive testing and validation
3. Performance optimization
4. Prepare for scientific applications

---

## 10. References and Resources

### Key Files
- `model/vegetation/soil.py` - Current empirical soil model
- `model/microbiome/DEMENTpy/src/dementpy.py` - DEMENTpy main
- `model/microbiome/DEMENTpy/src/grid.py` - Grid class
- `model/vegetation/gappy.py` - GAPPY main model

### Documentation
- GAPPY original papers (Shugart et al.)
- DEMENTpy publication (Wang et al.)
- CLAUDE.md - Project overview
- This document - Integration plan

---

## Appendices

### Appendix A: Code Snippets

**Example adapter initialization:**
```python
# In model/vegetation/gappy.py
from integration.dement_adapter import DEMENTpyAdapter

class GAPPYModel:
    def __init__(self, use_dement=False):
        # ... existing init ...
        self.use_dement = use_dement

        if use_dement:
            dement_config = {
                'pulse': 1,
                'end_time': 365,
                'interval': 30,
                'dispersal': 0
            }
            self.dement_adapter = DEMENTpyAdapter(dement_config)
```

### Appendix B: Unit Test Example

```python
# tests/test_integration.py
import pytest
from integration.dement_adapter import DEMENTpyAdapter

def test_litter_conversion():
    adapter = DEMENTpyAdapter(config)

    # GAPPY litter (tc/ha/day)
    litter_c1 = 0.01  # 10 kg C/ha/day

    # Convert to DEMENTpy (g/m²)
    substrate = adapter._litter_to_substrates(litter_c1, 0, 0, 0)

    # Should be: 0.01 tc/ha * 0.1 = 0.001 g/m² = 1 mg/m²
    assert abs(substrate['C'] - 1.0) < 1e-6
```

---

**Document End**
