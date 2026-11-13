Model and conventions
=====================

Geometry
--------

The sample is a height map :math:`z = h(x, y)` on a regular grid
(:class:`~bse_recon.SurfaceGrid`). Detectors are triangle meshes
(:class:`~bse_recon.Detector`). The beam travels in a straight line from
:attr:`ElectronSource.position <bse_recon.ElectronSource>` to each pixel. All
lengths must use the same unit.

Forward model
-------------

For a pixel at :math:`\mathbf p` with unit surface normal :math:`\hat{\mathbf n}`,
the signal collected by one detector is

.. math::

   S(\mathbf p) = \sum_{\text{triangles } t}
       g\bigl(\cos\theta_\text{in},\ \cos\theta_{\text{out},t},\ \cos\theta_{\text{spec},t}\bigr)\;
       \Omega_t(\mathbf p)

- :math:`\Omega_t` is the exact solid angle of triangle :math:`t`
  (Van Oosterom & Strackee, 1983).
- :math:`\theta_\text{in}` is the angle between :math:`\hat{\mathbf n}` and
  the direction back to the source.
- :math:`\theta_{\text{out},t}` is the angle between :math:`\hat{\mathbf n}`
  and the direction to the triangle centroid.
- :math:`\theta_{\text{spec},t}` is the angle between the mirror reflection
  of the beam and the direction to the triangle centroid.
- :math:`g` is the **scattering model**. The default,
  :func:`~bse_recon.lambertian`, is :math:`g = \eta_0 \max(\cos\theta_\text{out}, 0)/\pi`.
  Any callable with the same signature can be passed as ``model=``.

With ``normalize=True`` (the default), each image is divided by the image of
a flat surface on the same grid. This is the simulated equivalent of
flat-field correction, and it removes the shading caused by the beam moving
relative to the detectors.

Not modelled: shadowing by the surface, multiple scattering, and detector
efficiency or energy response.

Backscatter coefficient
~~~~~~~~~~~~~~~~~~~~~~~

For a pure element, :func:`~bse_recon.eta0_element` gives the
normal-incidence coefficient as a cubic in atomic number :math:`Z`. For an
alloy, :class:`~bse_recon.Material` takes the mass-weighted mean. The
dependence on incidence angle follows Arnal et al. (1969):

.. math::

   \eta(\theta) = B\,(\eta_0 / B)^{\cos\theta}, \qquad B = 0.89

Sign conventions
----------------

=====================  ============================================================================
Quantity               Convention
=====================  ============================================================================
Surface normal         :math:`(-\partial h/\partial x,\ -\partial h/\partial y,\ 1)`, normalised.
                       A surface rising towards :math:`+x` faces :math:`-x`.
Detector order         ``[left (-x), top (+y), right (+x), bottom (-y)]``
Contrast               :math:`M_x = (L-R)/(L+R)`, :math:`M_y = (B-T)/(B+T)`
=====================  ============================================================================

Slope inversion
---------------

Treat each detector as a point at :math:`(\pm d_w, 0, d_h)` or
:math:`(0, \pm d_w, d_h)`, and assume Lambertian emission. The contrasts at
beam position :math:`(x, y)` are then

.. math::

   M_x = \frac{g_x d_w}{g_x x + g_y y + d_h}, \qquad
   M_y = \frac{g_y d_w}{g_x x + g_y y + d_h},

which :func:`~bse_recon.slopes_from_contrast` inverts exactly:

.. math::

   g_x = \frac{M_x d_h}{d_w - M_x x - M_y y}, \qquad
   g_y = \frac{M_y d_h}{d_w - M_x x - M_y y}.

Two effects make this approximate in practice:

- **Surface relief.** A relief of :math:`\Delta z` causes a relative slope
  error of about :math:`\Delta z / d_h`.
- **Finite detector size.** Large detectors close to the sample break the
  point-detector assumption. See :doc:`auto_examples/05_detector_design_sweep`.

If the detector set is rotated about the beam axis by an angle ``azimuth``,
pass the same angle to :func:`~bse_recon.slopes_from_contrast`.

Height reconstruction
---------------------

:func:`~bse_recon.integrate_gradients` minimises
:math:`\lVert D_y Z - G_y\rVert^2 + \lVert Z D_x^\top - G_x\rVert^2`
(optionally :math:`+\lambda^2\lVert Z\rVert^2`). Here :math:`D_x` and
:math:`D_y` are second-order difference matrices. The minimiser is found
exactly by solving a Sylvester equation (Harker & O'Leary, 2008). Slopes
cannot determine the absolute height, so the result has zero mean.

Using measured images
---------------------

Pass real images in the same physical order (left, top, right, bottom in
the sample coordinates), with the same ``offset`` (:math:`d_w`) and
``height`` (:math:`d_h`). If one detector pair is mislabelled, the sign of
that slope component flips.

References
----------

- Arnal, F.; Verdier, P.; Vincensini, P.-D. *C. R. Acad. Sci. Paris* **268**, 1526 (1969).
- Harker, M.; O'Leary, P. *Least squares surface reconstruction from measured gradient fields.* CVPR (2008).
- Renner, J. et al. *A Ray Tracing Model for Electron Optical Imaging in Electron Beam Powder Bed Fusion.*
  J. Manuf. Mater. Process. **7**, 87 (2023). https://doi.org/10.3390/jmmp7030087
- Van Oosterom, A.; Strackee, J. *The Solid Angle of a Plane Triangle.* IEEE Trans. Biomed. Eng. **30**, 125 (1983).
