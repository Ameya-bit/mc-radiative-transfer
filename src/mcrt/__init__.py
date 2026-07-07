"""Monte Carlo radiative transfer through a plane-parallel neutron star atmosphere."""

from .beaming import (
    beaming_lookup,
    extract_intensity,
    fit_limb_darkening,
    fit_limb_darkening_slope,
)
from .convergence import (
    find_knee,
    loglog_slope,
    n_for_target_error,
    statistical_error,
)
from .monte_carlo import Photon, Simulation
from .pulse import (
    PulseProfile,
    analytic_isotropic_pf,
    bend,
    compute_profile,
    cos_psi,
    point_spot_flux,
    pulsed_fraction,
    visibility_threshold,
)
from .theory import chandrasekhar_h, eddington_limb_darkening
from .utils import (
    get_random_direction,
    rotate_vector,
    sample_step_size,
    sample_thomson_angle,
)

__all__ = [
    "Photon",
    "PulseProfile",
    "Simulation",
    "analytic_isotropic_pf",
    "beaming_lookup",
    "bend",
    "chandrasekhar_h",
    "compute_profile",
    "cos_psi",
    "eddington_limb_darkening",
    "extract_intensity",
    "find_knee",
    "fit_limb_darkening",
    "fit_limb_darkening_slope",
    "get_random_direction",
    "loglog_slope",
    "n_for_target_error",
    "point_spot_flux",
    "pulsed_fraction",
    "rotate_vector",
    "sample_step_size",
    "sample_thomson_angle",
    "statistical_error",
    "visibility_threshold",
]
