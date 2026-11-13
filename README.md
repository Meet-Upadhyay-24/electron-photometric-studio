# bse-recon

**Simulate backscattered-electron (BSE) photometric stereo: from a surface and a detector layout to synthetic detector images, slope maps and reconstructed heights.**

In electron-beam powder-bed fusion (and SEM), a focused beam scans the surface and BSE detectors placed around the beam record how many electrons reach them. Electrons leave a tilted surface preferentially in the direction it faces, so the difference between opposed detectors encodes the local slope. Integrating the slopes recovers the surface height. `bse-recon` lets you simulate this loop for your own detector geometry, material and surface, so you can design a layout or estimate measurement error before building hardware.

```
surface ─► simulate_images ─► contrast ─► slopes_from_contrast ─► integrate_gradients ─► height map
           (per detector)      Mx, My        ∂h/∂x, ∂h/∂y            (Sylvester LSQ)
```

## Install

```bash
git clone <repo-url>
cd bse_recon
pip install -e .            # numpy, scipy, matplotlib
pip install -e ".[dev]"     # + pytest, ruff
```

## Quick start

```python
import bse_recon as br

# geometry — all lengths in one unit (mm here)
detectors = br.four_quadrant_detectors(diameter=80, tilt=75, offset=132, height=87)
source    = br.ElectronSource(position=[0, 0, 1600])        # beam pivot point
surface   = br.SurfaceGrid.from_extent(150, 150, (-20, 20), (-20, 20))
surface   = surface.with_heights(br.sine_surface(surface, amplitude=0.25, wavenumber=2))

# forward model → four images [left, top, right, bottom]
images = br.simulate_images(surface, detectors, source)

# inversion
gx, gy = br.estimate_slopes(images, surface.X, surface.Y, offset=132, height=87)
z      = br.integrate_gradients(gx, gy, surface.dx, surface.dy)

print(br.height_errors(surface.heights, z))   # {'rmse': ..., 'mae': ..., 'max': ...}
```

## Examples

| Script | What it shows |
|---|---|
| [`01_detector_geometry.py`](examples/01_detector_geometry.py) | Shipped, generated and `.inp` detector meshes; placing them; the 4-detector layout |
| [`02_simulate_and_reconstruct.py`](examples/02_simulate_and_reconstruct.py) | The full loop with images, slope and height plots |
| [`03_beam_position_study.py`](examples/03_beam_position_study.py) | Error vs. scan-window position (the publication's "origins" study) |
| [`04_materials_and_scattering_models.py`](examples/04_materials_and_scattering_models.py) | Custom alloys and user-defined scattering models |
| [`05_detector_design_sweep.py`](examples/05_detector_design_sweep.py) | Sweep tilt / height / diameter, write a CSV of errors |
| [`06_your_own_surface.py`](examples/06_your_own_surface.py) | Simulate images of your own measured height map (`.npy`/`.csv`) |

Each script writes its figures to `examples/output/`.

## Model and conventions

**Geometry.** The sample is a height map `z = h(x, y)` on a regular grid (`SurfaceGrid`). Detectors are triangle meshes (`Detector`). The beam travels in a straight line from `ElectronSource.position` to each pixel.

**Forward model.** For each pixel `p` with unit normal `n`:

$$S(p) = \sum_{\text{triangles}} g(\cos\theta_\text{in}, \cos\theta_\text{out}, \cos\theta_\text{spec})\;\Omega_\text{tri}(p)$$

- `Ω` is the exact solid angle of each triangle (Van Oosterom & Strackee).
- `g` is the scattering model. The default, `lambertian`, is `η0 · max(cos θ_out, 0) / π`. Pass any callable as `model=` to change it (see `bse_recon.scattering`).
- With `normalize=True` (the default), each image is divided by the image of a flat surface on the same grid. This is the simulated equivalent of flat-field correction.
- Shadowing, multiple scattering and detector efficiency are not modelled.

**Sign conventions.**

| Quantity | Convention |
|---|---|
| Surface normal | `(−∂h/∂x, −∂h/∂y, 1)`, normalised. A surface rising towards +x faces −x. |
| Detector order | `[left (−x), top (+y), right (+x), bottom (−y)]`, as returned by `four_quadrant_detectors` |
| Contrast | `Mx = (L − R)/(L + R)`, `My = (B − T)/(B + T)` |
| Slopes | `gx = Mx·h / (w − Mx·x − My·y)`, `gy = My·h / (w − Mx·x − My·y)` with `w = offset`, `h = height` |

The slope formula is exact for small detectors when the surface relief is much smaller than `height`. Relief `Δz` adds a relative slope error of about `Δz / height`. Finite detector size adds the error that example 05 quantifies. If the detector set is rotated about the beam axis, pass the same `azimuth` to `slopes_from_contrast`.

**Applying it to measured images.** Pass real images in the same physical order (left, top, right, bottom in your sample coordinates) and use the same `offset` and `height`. A mislabelled pair flips the sign of that slope component.

## Validation

`pytest` runs 51 tests:

- **Against the original research code.** The original code and the notebook that produced the publication figures are reproduced for 8 scenarios: beam-position, detector-asymmetry and detector-diameter studies. When fed the original's inputs, the package reproduces the original detector images to float32 precision and the original Sylvester reconstruction to 1×10⁻⁷. With its own conventions it reaches the same reconstruction error within 5%. At the published 300×300 resolution the MAE values match within 2%: 0.0248 vs 0.0245, 0.182 vs 0.185 and 0.343 vs 0.344 mm. `tests/data/make_*_reference.py` documents how the reference data were made.
- **Analytic checks.** A closed triangulated surface subtends 4π, an on-axis disk matches its closed-form solid angle, mirror symmetry holds, and slopes of tilted planes are recovered at several azimuths. The closed-form inversion is exact for point detectors, and the integrator is exact for exact gradients.

```bash
pytest
```

## API overview

| Module | Contents |
|---|---|
| `bse_recon.surface` | `SurfaceGrid`, `sine_surface`, `gaussian_bump`, `random_surface`, `edge_taper` |
| `bse_recon.detector` | `Detector` (`from_inp`, `from_mesh`, `disk`, `rectangle`, `default_disk`, transforms), `four_quadrant_detectors` |
| `bse_recon.materials` | `Material`, `TI64`, `eta0_element` |
| `bse_recon.scattering` | `lambertian`, `lambertian_incidence`, and the model interface |
| `bse_recon.simulate` | `simulate_images`, `simulate_signal`, `detector_signal_arrays`, `triangle_solid_angles` |
| `bse_recon.gradients` | `contrast`, `slopes_from_contrast`, `estimate_slopes` |
| `bse_recon.reconstruction` | `integrate_gradients`, `remove_background`, `height_errors` |
| `bse_recon.plotting` | `plot_setup`, `plot_detector_2d`, `plot_images`, `plot_slopes`, `plot_reconstruction` |

## Citing

If you use this package, please cite it using [`CITATION.cff`](CITATION.cff). The publication reference will be added there.

## License

MIT — see [LICENSE](LICENSE).

## Acknowledgments and references

- Renner, J.; Grund, J.; Markl, M.; Körner, C. *A Ray Tracing Model for Electron Optical Imaging in Electron Beam Powder Bed Fusion.* J. Manuf. Mater. Process. 2023, 7, 87. <https://doi.org/10.3390/jmmp7030087>
- Van Oosterom, A.; Strackee, J. *The Solid Angle of a Plane Triangle.* IEEE Trans. Biomed. Eng. 1983, BME-30, 125–126.
- Arnal, F.; Verdier, P.; Vincensini, P.-D. *Coefficient de rétrodiffusion dans le cas d'électrons monocinétiques arrivant sur la cible sous une incidence oblique.* C. R. Acad. Sci. Paris 1969, 268, 1526–1529.
- Harker, M.; O'Leary, P. *Least squares surface reconstruction from measured gradient fields.* CVPR 2008.
