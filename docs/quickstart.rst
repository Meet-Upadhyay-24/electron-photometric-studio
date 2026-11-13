Quick start
===========

.. code-block:: python

   import bse_recon as br

   # geometry: all lengths in one unit (mm here)
   detectors = br.four_quadrant_detectors(diameter=80, tilt=75, offset=132, height=87)
   source    = br.ElectronSource(position=[0, 0, 1600])        # beam pivot point
   surface   = br.SurfaceGrid.from_extent(150, 150, (-20, 20), (-20, 20))
   surface   = surface.with_heights(br.sine_surface(surface, amplitude=0.25, wavenumber=2))

   # forward model: four images [left, top, right, bottom]
   images = br.simulate_images(surface, detectors, source)

   # inversion
   gx, gy = br.estimate_slopes(images, surface.X, surface.Y, offset=132, height=87)
   z      = br.integrate_gradients(gx, gy, surface.dx, surface.dy)

   print(br.height_errors(surface.heights, z))   # {'rmse': ..., 'mae': ..., 'max': ...}

The steps:

1. **Geometry.** A :class:`~bse_recon.SurfaceGrid` holds the height map. A
   :class:`~bse_recon.Detector` is a triangle mesh, and
   :func:`~bse_recon.four_quadrant_detectors` builds the usual symmetric set of
   four. The :class:`~bse_recon.ElectronSource` is the point the beam pivots
   about.
2. **Forward model.** :func:`~bse_recon.simulate_images` computes the signal
   each detector collects as the beam scans every pixel.
3. **Inversion.** :func:`~bse_recon.estimate_slopes` converts the four images
   into slopes, and :func:`~bse_recon.integrate_gradients` turns the slopes
   into heights.

The :doc:`auto_examples/index` go through each step with figures. Every
example can be downloaded as a Jupyter notebook.
