"""Monte Carlo radiative transfer through a plane-parallel neutron star atmosphere."""

from .beaming import extract_intensity, fit_limb_darkening_slope
from .convergence import (
    find_knee,
    loglog_slope,
    n_for_target_error,
    statistical_error,
)
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
    "find_knee",
    "fit_limb_darkening_slope",
    "get_random_direction",
    "loglog_slope",
    "n_for_target_error",
    "rotate_vector",
    "sample_step_size",
    "sample_thomson_angle",
    "statistical_error",
]
