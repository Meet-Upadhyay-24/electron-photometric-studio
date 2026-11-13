bse-recon
=========

**Simulate backscattered-electron (BSE) photometric stereo**: from a surface
and a detector layout to synthetic detector images, slope maps and
reconstructed heights.

In electron-beam powder-bed fusion and in the SEM, a focused beam scans the
surface while BSE detectors around the beam record how many electrons reach
them. A tilted surface sends more electrons in the direction it faces, so the
difference between opposed detectors encodes the local slope. Integrating
the slopes recovers the height. ``bse-recon`` simulates this loop for your own
detector geometry, material and surface, so you can design a detector layout
or estimate measurement error before building hardware.

.. code-block:: text

   surface ─► simulate_images ─► contrast ─► slopes_from_contrast ─► integrate_gradients ─► height map

- :doc:`installation`: install the package
- :doc:`quickstart`: the full loop in ten lines
- :doc:`auto_examples/index`: worked examples, each downloadable as a script or Jupyter notebook
- :doc:`model`: the physics, assumptions and sign conventions
- :doc:`api`: every public function and class

.. toctree::
   :maxdepth: 2
   :hidden:

   installation
   quickstart
   auto_examples/index
   model
   validation
   api
   changelog
