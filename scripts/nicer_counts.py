"""Adopted NICER count levels for PSR J0030+0451 and PSR J0740+6620.

Source of the numbers
---------------------
Miller et al. (2019, ApJL 887, L24) Fig. 10 and Miller et al. (2021, ApJL 918,
L28) Fig. 6 plot the bolometric counts per phase bin (black), the fitted
unmodulated background (dashed blue) and the full model (red). Both panels
were digitised from 400-dpi renders of the published PDFs: the red model
line was traced (it overlays the data bin-for-bin), the y axis was verified
against the printed tick marks, and the two plotted cycles were folded. The
per-bin profiles live in ``data/nicer_counts_j0030.csv`` and
``data/nicer_counts_j0740.csv`` (columns: bin, phase_center, total_counts,
spot_counts, background_counts).

"Spot counts" are the *modulated* counts: Miller's free, phase-independent
background per energy channel absorbs all unmodulated emission, so the real
pulsar contributes more than the spot profile shows (Bogdanov et al. 2019,
Table 2, quote a J0030 source rate of 314 counts/ks against the ~155 counts/ks
of modulated signal here).

Cross-checks
------------
* J0030: one-cycle total 1.34 M counts over 1.937 Ms = 0.69 counts/s, against
  the ~0.7 counts/s quoted in Miller et al. (2019) §2.2. PF(total) = 0.20.
* J0740: one-cycle total 549 k counts against the 521,004 events reported in
  Miller et al. (2021) §2.1 after filtering — the figure sits ~5 % above the
  text; the calibration is exact, so this is internal to that paper.
* J0740 with the NICER 3C50 background model (Salmi et al. 2022, ApJ 941, 150,
  Fig. 4): total 0.368–0.398 counts/s, background 0.354 counts/s, so
  B/<S> ≈ 12. Used as the high-background bound.

This module holds constants only. It is read by the revision scripts
(background ΔPF and the counts-level χ² test); the submitted anchor scripts
do not depend on it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parents[1] / "data"

# Salmi et al. (2022) 3C50 data set, 0.3–1.5 keV, 1.556 Ms: background
# 0.354 counts/s against a mean net source rate of ≈0.03 counts/s.
J0740_3C50_BACKGROUND_TO_MEAN_SPOT = 0.354 / 0.03


@dataclass(frozen=True)
class NicerCounts:
    """Phase-resolved count levels adopted from one NICER analysis figure."""

    name: str
    provenance: str
    exposure_s: float
    band_kev: tuple[float, float]
    phase: np.ndarray  # bin centres, cycles
    spot_counts: np.ndarray  # modulated counts per bin
    background_per_bin: float  # fitted unmodulated background, same for every bin

    @property
    def n_bins(self) -> int:
        return len(self.spot_counts)

    @property
    def total_counts(self) -> np.ndarray:
        return self.spot_counts + self.background_per_bin

    @property
    def background_to_mean_spot(self) -> float:
        return float(self.background_per_bin / self.spot_counts.mean())


def pulsed_fraction(counts: np.ndarray) -> float:
    """(max − min)/(max + min) of a phase-binned profile."""
    lo, hi = float(np.min(counts)), float(np.max(counts))
    return (hi - lo) / (hi + lo)


def _read_only(a: np.ndarray) -> np.ndarray:
    out = np.array(a, dtype=float)
    out.flags.writeable = False
    return out


def _load(csv_name: str, **fields) -> NicerCounts:
    table = np.genfromtxt(DATA_DIR / csv_name, delimiter=",", names=True)
    background = float(table["background_counts"][0])
    if not np.all(table["background_counts"] == background):
        raise ValueError(f"{csv_name}: background column is not constant")
    return NicerCounts(
        phase=_read_only(table["phase_center"]),
        spot_counts=_read_only(table["spot_counts"]),
        background_per_bin=background,
        **fields,
    )


SOURCES: dict[str, NicerCounts] = {
    "J0030": _load(
        "nicer_counts_j0030.csv",
        name="PSR J0030+0451",
        provenance="Miller et al. (2019) Fig. 10, digitised",
        exposure_s=1_936_864.0,
        band_kev=(0.4, 3.0),
    ),
    "J0740": _load(
        "nicer_counts_j0740.csv",
        name="PSR J0740+6620",
        provenance="Miller et al. (2021) Fig. 6, digitised",
        exposure_s=1_602_680.0,
        band_kev=(0.3, 1.23),
    ),
}


if __name__ == "__main__":
    for key, s in SOURCES.items():
        print(f"{key}: {s.provenance}; {s.n_bins} bins, {s.band_kev} keV, {s.exposure_s:.0f} s")
        print(f"  background/bin {s.background_per_bin:.0f}, mean spot/bin {s.spot_counts.mean():.0f}, "
              f"B/<S> = {s.background_to_mean_spot:.2f}")
        print(f"  PF(total) = {pulsed_fraction(s.total_counts):.3f}, "
              f"PF(spot) = {pulsed_fraction(s.spot_counts):.3f}")
