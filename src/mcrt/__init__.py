"""Monte Carlo radiative transfer through a plane-parallel neutron star atmosphere."""

from .beaming import extract_intensity, fit_limb_darkening_slope
from .monte_carlo import Photon, Simulation
from .theory import chandrasekhar_h, eddington_limb_darkening
from .utils import (
    get_random_direction,
    rotate_vector,
    sample_step_size,
    sample_thomson_angle,
)

__all__ = [
    "Photon",
    "Simulation",
    "chandrasekhar_h",
    "eddington_limb_darkening",
    "extract_intensity",
    "fit_limb_darkening_slope",
    "get_random_direction",
    "rotate_vector",
    "sample_step_size",
    "sample_thomson_angle",
]
