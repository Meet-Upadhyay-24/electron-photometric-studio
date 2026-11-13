Validation
==========

The test suite (``pytest``) checks the package in two ways.

Against the original research code
----------------------------------

The original code, and the notebook that produced the publication figures,
were rerun for 8 scenarios: the beam-position study, the detector-asymmetry
study and the detector-diameter study. Their outputs are stored in
``tests/data``, and the scripts that generated them sit alongside.

- Given the original's inputs, the package reproduces the original detector
  images to float32 precision and the Sylvester reconstruction to
  :math:`10^{-7}`.
- With its own conventions, which are listed in the :doc:`changelog`, the
  package reaches the same reconstruction error to within 5 %.

At the published 300 × 300 resolution, the beam-position study gives these
mean absolute errors (offset removed):

=================  ===============  =========
Scan centre (mm)   Original code    Package
=================  ===============  =========
(0, 0)             0.0245 mm        0.0248 mm
(40, 40)           0.185 mm         0.182 mm
(80, 80)           0.344 mm         0.343 mm
=================  ===============  =========

Against analytic results
------------------------

- A closed triangulated surface subtends :math:`4\pi`.
- A disk seen on its axis matches the closed-form solid angle.
- A flat surface gives contrasts that are antisymmetric about the beam axis.
- Slopes of tilted planes are recovered at several detector azimuths.
- The closed-form inversion is exact for point detectors.
- Integrating exact gradients reproduces the surface.
