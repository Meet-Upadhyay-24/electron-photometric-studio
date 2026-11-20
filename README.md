# bse-recon

**Simulate backscattered-electron (BSE) photometric stereo: from a surface and a detector layout to synthetic detector images, slope maps and reconstructed heights.**

In electron-beam powder-bed fusion (and SEM), a focused beam scans the surface and BSE detectors placed around the beam record how electron intensity. Electrons leave a tilted surface preferentially in the direction it faces, so the difference between opposed detectors encodes the local slope. Integrating the slopes recovers the surface height. `bse-recon` lets you simulate this loop for your own detector geometry, material and surface, so you can design a layout or estimate measurement error before building hardware.

```
surface ─► simulate_images ─► contrast ─► slopes_from_contrast ─► integrate_gradients ─► height map
           (per detector)      Mx, My        ∂h/∂x, ∂h/∂y            (Sylvester LSQ)
```

**Documentation:** <https://meet-upadhyay-24.github.io/electron-photometric-studio/>. It covers the model, the API reference and a worked example gallery. Each example can be downloaded as a Jupyter notebook.

## Install

```bash
git clone https://github.com/Meet-Upadhyay-24/electron-photometric-studio.git
cd electron-photometric-studio
pip install -e .            # numpy, scipy, matplotlib
pip install -e ".[docs]"    # + Sphinx, to build the documentation
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

The [example gallery](https://meet-upadhyay-24.github.io/electron-photometric-studio/auto_examples/index.html) shows each example rendered with its figures, plus `.py` and Jupyter notebook downloads.

| Script | What it shows |
|---|---|
| [`01_detector_geometry.py`](examples/01_detector_geometry.py) | Shipped, generated and `.inp` detector meshes; placing them; the 4-detector layout |
| [`02_simulate_and_reconstruct.py`](examples/02_simulate_and_reconstruct.py) | The full loop with images, slope and height plots |
| [`03_beam_position_study.py`](examples/03_beam_position_study.py) | Error vs. scan-window position |
| [`04_materials_and_scattering_models.py`](examples/04_materials_and_scattering_models.py) | Custom alloys and user-defined scattering models |
| [`05_detector_design_sweep.py`](examples/05_detector_design_sweep.py) | Sweep tilt / height / diameter, write a CSV of errors |
| [`06_your_own_surface.py`](examples/06_your_own_surface.py) | Simulate images of your own measured height map (`.npy`/`.csv`) |

Run any script directly (`python examples/02_simulate_and_reconstruct.py`) to open its figures. Examples 05 and 06 also save data to an `output/` folder.

## Model and conventions

**Geometry.** The sample is a height map `z = h(x, y)` on a regular grid (`SurfaceGrid`). Detectors are triangle meshes (`Detector`). The beam travels in a straight line from `ElectronSource.position` to each grid point.

**Forward model.** For each pixel `p` with unit normal `n`:

$$S(p) = \sum_{\text{triangles}} g(\cos\theta_\text{in}, \cos\theta_\text{out}, \cos\theta_\text{spec})\;\Omega_\text{tri}(p)$$

- `Ω` is the exact solid angle of each triangle (Van Oosterom & Strackee).
- `g` is the scattering model. The default, `lambertian`, is `η0 · max(cos θ_out, 0) / π`. Pass any callable as `model=` to change it (see `bse_recon.scattering`).
- With `normalize=True` (the default), each image is divided by the image of a flat surface on the same grid. This is the simulated equivalent of flat-field correction.
- Shadowing, multiple scattering and detector efficiency are not modelled.

**Sign conventions.**

| Quantity | Convention |
|---|---|
| Surface normal | `(−∂h/∂x, −∂h/∂y, 1)`, normalised. |
| Detector order | `[left (−x), top (+y), right (+x), bottom (−y)]`, as returned by `four_quadrant_detectors` |
| Contrast | `Mx = (L − R)/(L + R)`, `My = (B − T)/(B + T)` |
| Slopes | `gx = Mx·h / (w − Mx·x − My·y)`, `gy = My·h / (w − Mx·x − My·y)` with `w = offset`, `h = height` |

The slope formula is exact for small detectors when the surface relief is much smaller than `height`. Relief `Δz` adds a relative slope error of about `Δz / height`.

**Applying it to measured images.** Pass real images in the same physical order (left, top, right, bottom in your sample coordinates) and use the same `offset` and `height`. A mislabelled pair flips the sign of that slope component.


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

If you use this package, please cite the accompanying article:

> Upadhyay, M.; McDonald, T.; Horn, T.; Britton, T. B.; Sinclair, C.; Maijer, D. *Enhancing backscatter electron based surface height reconstructions in electron beam powder bed fusion and exploring correlations to porosity.* Journal of Manufacturing Processes 2026, 172, 54–69. <https://doi.org/10.1016/j.jmapro.2026.05.034>

GitHub's **Cite this repository** button (from [`CITATION.cff`](CITATION.cff)) gives the same reference in BibTeX and APA formats.

## License

MIT — see [LICENSE](LICENSE).

## Acknowledgments and references

- Upadhyay, M.; McDonald, T.; Horn, T.; Britton, T. B.; Sinclair, C.; Maijer, D. *Enhancing backscatter electron based surface height reconstructions in electron beam powder bed fusion and exploring correlations to porosity.* J. Manuf. Process. 2026, 172, 54–69. <https://doi.org/10.1016/j.jmapro.2026.05.034>
- Renner, J.; Grund, J.; Markl, M.; Körner, C. *A Ray Tracing Model for Electron Optical Imaging in Electron Beam Powder Bed Fusion.* J. Manuf. Mater. Process. 2023, 7, 87. <https://doi.org/10.3390/jmmp7030087>
- Harker, M.; O'Leary, P. *Least squares surface reconstruction from measured gradient fields.* 2008 IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 1–7. <https://doi.org/10.1109/CVPR.2008.4587414>
