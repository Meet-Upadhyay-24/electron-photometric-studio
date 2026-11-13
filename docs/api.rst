API reference
=============

.. currentmodule:: bse_recon

Everything below can be imported from the top level: ``import bse_recon as br``.

Geometry
--------

.. autosummary::
   :toctree: generated
   :nosignatures:

   SurfaceGrid
   Detector
   four_quadrant_detectors
   ElectronSource

Synthetic surfaces
------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   sine_surface
   gaussian_bump
   random_surface
   edge_taper

Materials and scattering
------------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   Material
   eta0_element
   lambertian
   lambertian_incidence

The predefined material :data:`bse_recon.TI64` is Ti-6Al-4V.

Simulation
----------

.. autosummary::
   :toctree: generated
   :nosignatures:

   simulate_images
   simulate_signal
   detector_signal_arrays
   triangle_solid_angles

Inversion and reconstruction
----------------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   contrast
   slopes_from_contrast
   estimate_slopes
   integrate_gradients
   remove_background
   height_errors

Plotting
--------

.. currentmodule:: bse_recon.plotting

.. autosummary::
   :toctree: generated
   :nosignatures:

   plot_setup
   plot_detector_2d
   plot_images
   plot_slopes
   plot_reconstruction
   apply_science_style
