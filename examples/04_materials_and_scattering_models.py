"""
04 — Change the material and plug in your own scattering model.

A scattering model is any function ``g(cos_in, cos_out, cos_spec, material)``
(see ``bse_recon.scattering``). Here we compare:

* ``lambertian``            — cosine emission, constant yield (package default)
* ``lambertian_incidence``  — cosine emission, yield rising with beam incidence angle
* a user-defined model with a forward-peaked (specular) lobe

Note that ``lambertian_incidence`` gives the same slopes as ``lambertian``: the
incidence-dependent yield η(θ) scales all four detectors equally at each pixel,
so it cancels in the normalised contrast. A directional (specular) lobe does
not cancel and biases the estimate.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import bse_recon as br

OUT = Path(__file__).parent / "output"
OUT.mkdir(exist_ok=True)

steel = br.Material("316L", ((26, 0.68), (24, 0.17), (28, 0.12), (42, 0.03)))
print(f"η0  Ti-6Al-4V = {br.TI64.eta0:.3f}   316L = {steel.eta0:.3f}")


def lambert_plus_specular(cos_in, cos_out, cos_spec, material, weight=0.3, k=8):
    """Diffuse cosine lobe plus a cos^k lobe around the mirror direction."""
    diffuse = np.maximum(cos_out, 0) / np.pi
    lobe = (k + 1) / (2 * np.pi) * np.maximum(cos_spec, 0) ** k
    return material.eta0 * ((1 - weight) * diffuse + weight * lobe)


offset, height = 132.0, 87.0
detectors = br.four_quadrant_detectors(diameter=80, tilt=75, offset=offset, height=height)
source = br.ElectronSource([0, 0, 1600])
s = br.SurfaceGrid.from_extent(120, 120, (-20, 20), (-20, 20))
s = s.with_heights(br.gaussian_bump(s, amplitude=1.0, sigma=4.0))
gx_true, _ = s.gradients()

models = {"lambertian": br.lambertian,
          "lambertian_incidence": br.lambertian_incidence,
          "lambert + specular": lambert_plus_specular}

fig, ax = plt.subplots(figsize=(7, 4))
k = s.shape[0] // 2
ax.plot(s.x, gx_true[k], "k-", lw=2, label="true ∂h/∂x")
for name, model in models.items():
    images = br.simulate_images(s, detectors, source, material=steel, model=model)
    gx, gy = br.estimate_slopes(images, s.X, s.Y, offset, height)
    z = br.integrate_gradients(gx, gy, s.dx, s.dy)
    print(f"{name:22s} RMSE {br.height_errors(s.heights, z)['rmse']:.4f} mm")
    ax.plot(s.x, gx[k], "--", label=name)
ax.set(xlabel="x [mm]", ylabel="slope", title="Slope estimated under different scattering models (316L)")
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "04_models.png", dpi=150)
plt.show()
