"""
Integration layer for coupling GAPPY and DEMENTpy models.

This package provides the coupling interface between the forest gap model (GAPPY)
and the microbial decomposition model (DEMENTpy) to create the integrated GoGs
(Gap of Gaps) ecosystem model.
"""

from .dement_adapter import DEMENTpyAdapter
from .unit_conversions import UnitConverter
from .parameter_mapping import ParameterMapper

__all__ = ['DEMENTpyAdapter', 'UnitConverter', 'ParameterMapper']

__version__ = '0.1.0'
__author__ = 'GoGs Development Team'
