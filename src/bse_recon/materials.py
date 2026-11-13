"""Material models: backscatter coefficients for pure elements and alloys."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np

#: Saturation backscatter coefficient at grazing incidence (Arnal et al., 1969).
B_SATURATION = 0.89


def eta0_element(Z) -> np.ndarray:
    """
    Normal-incidence backscatter coefficient of a pure element.

    Empirical cubic fit in atomic number ``Z`` (Reuter, 1972), valid for
    beam energies of roughly 10-100 keV.
    """
    Z = np.asarray(Z, dtype=float)
    return -0.0254 + 0.016 * Z - 1.86e-4 * Z**2 + 8.3e-7 * Z**3


@dataclass(frozen=True)
class Material:
    """
    A material defined by its elemental composition.

    Parameters
    ----------
    name:
        Human-readable label.
    composition:
        Sequence of ``(atomic_number, mass_fraction)`` pairs; fractions must sum to 1.

    Examples
    --------
    >>> steel_316 = Material("316L", ((26, 0.68), (24, 0.17), (28, 0.12), (42, 0.03)))
    >>> round(steel_316.eta0, 3)
    0.282
    """

    name: str
    composition: Tuple[Tuple[int, float], ...]

    def __post_init__(self):
        object.__setattr__(self, "composition", tuple((int(z), float(f)) for z, f in self.composition))
        total = sum(f for _, f in self.composition)
        if not np.isclose(total, 1.0, atol=1e-6):
            raise ValueError(f"Mass fractions of {self.name!r} must sum to 1 (got {total:.6f}).")

    @property
    def eta0(self) -> float:
        """Normal-incidence backscatter coefficient (mass-fraction weighted mean)."""
        Z = np.array([z for z, _ in self.composition])
        f = np.array([f for _, f in self.composition])
        return float(np.sum(eta0_element(Z) * f))

    def eta(self, cos_incidence) -> np.ndarray:
        """
        Backscatter coefficient at an incidence angle θ from the surface normal.

        Uses ``η(θ) = B · (η0 / B) ** cos θ`` (Arnal et al., 1969), which gives
        ``η0`` at normal incidence and tends to ``B = 0.89`` at grazing incidence.

        Parameters
        ----------
        cos_incidence:
            ``cos θ``; values are clipped to ``[0, 1]``.
        """
        c = np.clip(np.asarray(cos_incidence, dtype=float), 0.0, 1.0)
        return B_SATURATION * (self.eta0 / B_SATURATION) ** c


#: Ti-6Al-4V (mass fractions).
TI64 = Material("Ti-6Al-4V", ((22, 0.90), (13, 0.06), (23, 0.04)))
