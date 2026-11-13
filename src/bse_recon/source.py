"""Electron beam source."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

E_CHARGE = 1.602176634e-19  # C


@dataclass
class ElectronSource:
    """
    Electron beam origin.

    Only ``position`` enters the geometric simulation: the beam is assumed to
    travel in a straight line from ``position`` to every surface node (i.e. the
    beam pivot point of a deflected beam). The remaining fields are optional
    metadata used to convert normalised signals into electron counts.

    Parameters
    ----------
    position:
        ``[x, y, z]`` of the beam pivot, in the same length unit as the surface.
    voltage:
        Accelerating voltage [V].
    current:
        Beam current [A].
    dwell_time:
        Dwell time per pixel [s].
    """

    position: Sequence[float]
    voltage: Optional[float] = None
    current: Optional[float] = None
    dwell_time: Optional[float] = None

    def __post_init__(self):
        self.position = np.asarray(self.position, dtype=float)
        if self.position.shape != (3,):
            raise ValueError("ElectronSource.position must be a 3-element vector.")

    @property
    def primary_electrons(self) -> float:
        """Primary electrons per pixel, ``current · dwell_time / e``."""
        if self.current is None or self.dwell_time is None:
            raise ValueError("Set current and dwell_time to compute the primary electron count.")
        return self.current * self.dwell_time / E_CHARGE
