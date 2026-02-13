"""
Unit conversion utilities for GAPPY-DEMENTpy coupling.

This module handles all unit conversions between the two models:
- GAPPY uses: tc/ha, tn/ha, cm, °C
- DEMENTpy uses: g/m², °C

Conversion factors:
- 1 tc/ha = 0.1 kg/m² = 100 g/m²
- 1 tn/ha = 0.1 kg/m² = 100 g/m²
"""

import numpy as np
from typing import Dict, Tuple


class UnitConverter:
    """Handles all unit conversions between GAPPY and DEMENTpy."""

    # Conversion factors
    TC_HA_TO_G_M2 = 100.0      # 1 tc/ha = 100 g/m²
    TN_HA_TO_G_M2 = 100.0      # 1 tn/ha = 100 g/m²
    G_M2_TO_TC_HA = 0.01       # 1 g/m² = 0.01 tc/ha
    G_M2_TO_TN_HA = 0.01       # 1 g/m² = 0.01 tn/ha

    # Time conversions
    DAY_TO_YEAR = 365.0
    YEAR_TO_DAY = 1.0 / 365.0

    @staticmethod
    def litter_to_substrate(litter_c_tc_ha: float, litter_n_tn_ha: float) -> Dict[str, float]:
        """
        Convert GAPPY litter inputs to DEMENTpy substrate pools.

        Args:
            litter_c_tc_ha: Litter carbon in tc/ha/day
            litter_n_tn_ha: Litter nitrogen in tn/ha/day

        Returns:
            Dictionary with 'C' and 'N' in g/m²/day
        """
        substrate_c = litter_c_tc_ha * UnitConverter.TC_HA_TO_G_M2
        substrate_n = litter_n_tn_ha * UnitConverter.TN_HA_TO_G_M2

        return {
            'C': substrate_c,
            'N': substrate_n,
            'C_N_ratio': substrate_c / substrate_n if substrate_n > 0 else 0.0
        }

    @staticmethod
    def annual_to_daily_litter(annual_c: float, annual_n: float) -> Tuple[float, float]:
        """
        Convert annual litter production to daily rates.

        Args:
            annual_c: Annual litter C in tc/ha/year
            annual_n: Annual litter N in tn/ha/year

        Returns:
            Tuple of (daily_c, daily_n) in tc/ha/day
        """
        daily_c = annual_c * UnitConverter.YEAR_TO_DAY
        daily_n = annual_n * UnitConverter.YEAR_TO_DAY

        return daily_c, daily_n

    @staticmethod
    def dement_to_gappy_n(dement_n_g_m2: float) -> float:
        """
        Convert DEMENTpy available N to GAPPY units.

        Args:
            dement_n_g_m2: Available N from DEMENTpy in g/m²

        Returns:
            Available N in tn/ha
        """
        return dement_n_g_m2 * UnitConverter.G_M2_TO_TN_HA

    @staticmethod
    def dement_to_gappy_resp(dement_resp_g_m2_day: float) -> float:
        """
        Convert DEMENTpy respiration to GAPPY units.

        Args:
            dement_resp_g_m2_day: CO2 respiration from DEMENTpy in g C/m²/day

        Returns:
            Respiration in tc/ha/day
        """
        return dement_resp_g_m2_day * UnitConverter.G_M2_TO_TC_HA

    @staticmethod
    def aggregate_plot_outputs(outputs_list: list, aggregation: str = 'mean') -> Dict[str, float]:
        """
        Aggregate outputs from multiple DEMENTpy grids (one per plot).

        Args:
            outputs_list: List of output dictionaries from each grid
            aggregation: Aggregation method ('mean', 'sum', 'median')

        Returns:
            Aggregated output dictionary
        """
        if not outputs_list:
            return {'avail_N': 0.0, 'C_resp': 0.0}

        # Stack arrays
        avail_n_array = np.array([out['avail_N'] for out in outputs_list])
        c_resp_array = np.array([out['C_resp'] for out in outputs_list])

        if aggregation == 'mean':
            agg_func = np.mean
        elif aggregation == 'sum':
            agg_func = np.sum
        elif aggregation == 'median':
            agg_func = np.median
        else:
            raise ValueError(f"Unknown aggregation method: {aggregation}")

        return {
            'avail_N': agg_func(avail_n_array),
            'C_resp': agg_func(c_resp_array),
            'avail_N_std': np.std(avail_n_array),
            'C_resp_std': np.std(c_resp_array)
        }

    @staticmethod
    def validate_mass_balance(inputs: Dict[str, float], outputs: Dict[str, float],
                            state_change: Dict[str, float], tolerance: float = 1e-6) -> bool:
        """
        Validate mass balance for C and N.

        Args:
            inputs: Dict with 'C_in' and 'N_in'
            outputs: Dict with 'C_out' and 'N_out'
            state_change: Dict with 'delta_C' and 'delta_N'
            tolerance: Acceptable relative error

        Returns:
            True if mass balance is satisfied within tolerance
        """
        # C balance: inputs = outputs + accumulation
        c_balance = inputs['C_in'] - outputs['C_out'] - state_change['delta_C']
        c_relative_error = abs(c_balance) / max(inputs['C_in'], 1e-10)

        # N balance: inputs = outputs + accumulation
        n_balance = inputs['N_in'] - outputs['N_out'] - state_change['delta_N']
        n_relative_error = abs(n_balance) / max(inputs['N_in'], 1e-10)

        c_ok = c_relative_error < tolerance
        n_ok = n_relative_error < tolerance

        if not (c_ok and n_ok):
            print(f"Mass balance warning:")
            print(f"  C error: {c_relative_error:.2e} (tolerance: {tolerance})")
            print(f"  N error: {n_relative_error:.2e} (tolerance: {tolerance})")

        return c_ok and n_ok


def print_conversion_table():
    """Print a reference table of unit conversions."""
    print("\n" + "="*60)
    print("GAPPY-DEMENTpy Unit Conversion Reference")
    print("="*60)
    print("\nCarbon:")
    print(f"  1 tc/ha = {UnitConverter.TC_HA_TO_G_M2} g/m²")
    print(f"  1 g/m² = {UnitConverter.G_M2_TO_TC_HA} tc/ha")
    print("\nNitrogen:")
    print(f"  1 tn/ha = {UnitConverter.TN_HA_TO_G_M2} g/m²")
    print(f"  1 g/m² = {UnitConverter.G_M2_TO_TN_HA} tn/ha")
    print("\nTime:")
    print(f"  1 year = {UnitConverter.DAY_TO_YEAR} days")
    print(f"  1 day = {UnitConverter.YEAR_TO_DAY} years")
    print("="*60 + "\n")


if __name__ == "__main__":
    # Demonstrate conversions
    print_conversion_table()

    # Example conversion
    litter_c = 0.01  # tc/ha/day
    litter_n = 0.001  # tn/ha/day

    substrate = UnitConverter.litter_to_substrate(litter_c, litter_n)
    print(f"Example: {litter_c} tc/ha/day → {substrate['C']:.4f} g/m²/day")
    print(f"Example: {litter_n} tn/ha/day → {substrate['N']:.4f} g/m²/day")
    print(f"C:N ratio: {substrate['C_N_ratio']:.2f}")
