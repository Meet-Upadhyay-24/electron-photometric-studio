# Changelog

## 0.1.0 — first public version

Restructured from the research script `main_arrays_parallel.py` and its helpers (`det_funcs.py`, `surf_reconstruction.py`).

### Verified against the research code
- Detector images reproduce the original to float32 precision, and the Sylvester reconstruction to 1×10⁻⁷.
- This was checked for the scenarios of the publication notebook (`paperfigures.ipynb`); see `tests/test_regression.py`.

### Changed relative to the research code
- **Surface normal** is the geometric outward normal `(−∂h/∂x, −∂h/∂y, 1)`. The research code used `(+∂h/∂x, +∂h/∂y, 1)` and compensated with `Mx = (R − L)/(R + L)` and `x_s = Mx·d_h/(d_w + Mx·x + My·y)`. The package uses `Mx = (L − R)/(L + R)` and `x_s = Mx·d_h/(d_w − Mx·x − My·y)`. Simulated slopes and heights agree with the original within a few per cent; only the sign of `Mx`/`My` is flipped.
- **Grid spacing** is `(max − min)/(n − 1)`, which is the actual `linspace` step. The research code used `(max − min)/n`, which over-stated simulated slopes by a factor `n/(n − 1)` (0.3 % at n = 300). That error cancelled in the reconstruction comparisons.
- **Detector vertices** are kept in float64 instead of float16. This changes images by about 2×10⁻⁵ relative.
- **Free-boundary Sylvester solve** adds a vanishing ridge so the singular system has a unique solution. Reconstructions are returned with zero mean.
  - The unregularised solve left an arbitrary constant offset. At origin (0,0) of the beam-position study, that offset raised the notebook's MAE from 0.025 to 0.035 mm.
- **Backscatter coefficient vs. incidence** follows Arnal et al. (1969), `η(θ) = B(η0/B)^cos θ`. The research code used `θ` in radians as the exponent. This never affected the published results, because the default Lambertian model uses a constant yield.
- **Absolute (non-normalised) signals** use `η0` as the yield. The research code used `η(90°)`, a factor of about 2.11 higher. Normalised images are unaffected.
- The unused specular term of the research `scattering_function` was dropped. User-defined models can add one (see example 04).
- The forward model processes detector triangles in chunks under a memory budget, so large grids no longer need `(ny, nx, n_triangles, 3)` arrays.
- Removed: `ProcessPoolExecutor` loop, `old_solid_angle`, `corrected_detector_image`, `plate_electrons`, `det_electrons`, the pylops dependency and hard-coded paths.
