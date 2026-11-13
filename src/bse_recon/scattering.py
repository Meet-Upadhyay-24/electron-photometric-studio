"""
Angular scattering models.

A scattering model is any callable with the signature::

    model(cos_in, cos_out, cos_spec, material) -> weight

returning the backscattered electrons per primary electron per steradian
emitted towards a detector element. All arguments broadcast against each other:

``cos_in``
    cosine of the angle between the surface normal and the direction from the
    surface point back to the beam source (1 at normal incidence).
``cos_out``
    cosine of the angle between the surface normal and the direction from the
    surface point to the detector element.
``cos_spec``
    cosine of the angle between the specular reflection of the incoming beam and
    the direction to the detector element.
``material``
    a :class:`bse_recon.Material`.

Negative weights are treated as zero by the simulator. Pass your own callable
to :func:`bse_recon.simulate_images` to try a different model.
"""
from __future__ import annotations

import numpy as np

from .materials import Material


def lambertian(cos_in, cos_out, cos_spec, material: Material) -> np.ndarray:
    """
    Cosine (Lambertian) emission with a constant total yield ``η0``.

    ``g = η0 · max(cos_out, 0) / π``.

    This is the model used for the published simulations: the detected signal
    depends on the surface tilt only through the emission direction.
    """
    return material.eta0 * np.maximum(cos_out, 0.0) / np.pi


def lambertian_incidence(cos_in, cos_out, cos_spec, material: Material) -> np.ndarray:
    """
    Cosine emission whose total yield follows the incidence-angle dependence ``η(θ)``.

    ``g = η(θ_in) · max(cos_out, 0) / π`` with ``η(θ)`` from :meth:`Material.eta`.
    Tilting the surface away from the beam increases the total backscatter yield.
    """
    return material.eta(cos_in) * np.maximum(cos_out, 0.0) / np.pi
