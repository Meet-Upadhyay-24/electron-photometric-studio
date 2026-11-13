"""
Materials and scattering models
===============================

How to change the material and plug in your own angular scattering model.
"""

import matplotlib.pyplot as plt
import numpy as np

import bse_recon as br

# %%
# Materials
# ---------
# A :class:`~bse_recon.Material` is a list of ``(atomic number, mass fraction)``
# pairs. Its normal-incidence backscatter coefficient ``η0`` is the
# mass-weighted mean of the elemental values.

steel = br.Material("316L", ((26, 0.68), (24, 0.17), (28, 0.12), (42, 0.03)))
print(f"eta0  Ti-6Al-4V = {br.TI64.eta0:.3f}   316L = {steel.eta0:.3f}")

# %%
# Scattering models
# -----------------
# A scattering model is any function ``g(cos_in, cos_out, cos_spec, material)``
# returning backscattered electrons per primary electron per steradian. Its
# arguments are the cosines of three angles: beam incidence, emission towards
# the detector element, and deviation from the mirror direction. The package
# provides:
#
# * :func:`~bse_recon.lambertian`: cosine emission with constant yield (default)
# * :func:`~bse_recon.lambertian_incidence`: cosine emission whose yield rises
#   with the beam incidence angle
#
# The model below adds a forward-peaked lobe around the mirror direction:


def lambert_plus_specular(cos_in, cos_out, cos_spec, material, weight=0.3, k=8):
    """Diffuse cosine lobe plus a cos^k lobe around the mirror direction."""
    diffuse = np.maximum(cos_out, 0) / np.pi
    lobe = (k + 1) / (2 * np.pi) * np.maximum(cos_spec, 0) ** k
    return material.eta0 * ((1 - weight) * diffuse + weight * lobe)


# %%
# Compare the models on a Gaussian hill
# -------------------------------------
# All three models are inverted with the same Lambertian slope formula, so any
# difference between their results shows how much model error biases the
# measurement.

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
ax.set(xlabel="x [mm]", ylabel="slope", title="Slope along the centre row (316L)")
ax.legend()

# %%
# ``lambertian_incidence`` gives exactly the same slopes as ``lambertian``.
# The incidence-dependent yield scales all four detectors equally at each
# pixel, so it cancels in the normalised contrast. A directional (specular)
# lobe does not cancel, and it biases the estimate.

plt.show()
