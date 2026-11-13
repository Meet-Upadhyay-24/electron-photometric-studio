Installation
============

``bse-recon`` needs Python 3.9 or newer, plus NumPy, SciPy and Matplotlib.

From source:

.. code-block:: bash

   git clone https://github.com/Meet-Upadhyay-24/electron-photometric-studio.git
   cd electron-photometric-studio
   pip install -e .

Optional extras:

.. code-block:: bash

   pip install -e ".[dev]"       # pytest and ruff, to run the test suite
   pip install -e ".[docs]"      # Sphinx, to build this documentation
   pip install -e ".[science]"   # SciencePlots figure style

Check the installation by running the test suite:

.. code-block:: bash

   pytest

Building the documentation
--------------------------

.. code-block:: bash

   pip install -e ".[docs]"
   sphinx-build -b html docs docs/_build/html

Every example is executed during the build, so a full build takes a few minutes.
