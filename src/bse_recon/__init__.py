"""
bse_recon — simulate backscattered-electron (BSE) photometric stereo.

A beam scans a surface; several BSE detectors record how many electrons reach
them. Because emission depends on the local surface tilt, differences between
opposed detectors encode the surface slopes, which can be integrated into a
height map. This package provides each step:

1. geometry   — :class:`SurfaceGrid`, :class:`Detector`, :func:`four_quadrant_detectors`, :class:`ElectronSource`
2. physics    — :class:`Material` and scattering models (:mod:`bse_recon.scattering`)
3. simulation — :func:`simulate_images`
4. inversion  — :func:`contrast`, :func:`slopes_from_contrast`, :func:`integrate_gradients`

Plotting helpers live in :mod:`bse_recon.plotting` (imported on demand).
"""
from .detector import Detector, four_quadrant_detectors
from .gradients import contrast, estimate_slopes, slopes_from_contrast
from .materials import TI64, Material, eta0_element
from .reconstruction import height_errors, integrate_gradients, remove_background
from .scattering import lambertian, lambertian_incidence
from .simulate import detector_signal_arrays, simulate_images, simulate_signal, triangle_solid_angles
from .source import ElectronSource
from .surface import SurfaceGrid, edge_taper, gaussian_bump, random_surface, sine_surface

__version__ = "0.1.0"

__all__ = [
    "Detector", "four_quadrant_detectors",
    "SurfaceGrid", "sine_surface", "gaussian_bump", "random_surface", "edge_taper",
    "ElectronSource",
    "Material", "TI64", "eta0_element",
    "lambertian", "lambertian_incidence",
    "simulate_images", "simulate_signal", "detector_signal_arrays", "triangle_solid_angles",
    "contrast", "slopes_from_contrast", "estimate_slopes",
    "integrate_gradients", "remove_background", "height_errors",
    "__version__",
]
