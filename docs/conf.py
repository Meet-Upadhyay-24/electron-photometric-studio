"""Sphinx configuration for bse-recon."""
import os
import sys

import matplotlib

matplotlib.use("Agg")
sys.path.insert(0, os.path.abspath("../src"))

import bse_recon  # noqa: E402

project = "bse-recon"
author = "Meet Upadhyay"
copyright = f"2026, {author}"
release = bse_recon.__version__
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_gallery.gen_gallery",
    "sphinx_copybutton",
    "myst_parser",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build"]
templates_path = ["_templates"]

# -- API reference ---------------------------------------------------------------
autosummary_generate = True
autodoc_default_options = {"members": True, "show-inheritance": False}
autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_numpy_docstring = True
napoleon_google_docstring = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
}

# -- Example gallery ---------------------------------------------------------------
sphinx_gallery_conf = {
    "examples_dirs": "../examples",
    "gallery_dirs": "auto_examples",
    "filename_pattern": r"[\\/]\d\d_",      # execute every numbered example (any OS path separator)
    "ignore_pattern": r"__init__\.py",
    "within_subsection_order": "FileNameSortKey",
    "download_all_examples": True,
    "doc_module": ("bse_recon",),
    "reference_url": {"bse_recon": None},
    "backreferences_dir": "gen_modules/backreferences",
    "remove_config_comments": True,
    "capture_repr": (),                     # show print() output, not "<Figure …>" reprs
    "matplotlib_animations": False,
    "image_scrapers": ("matplotlib",),
}

# -- HTML ------------------------------------------------------------------------
html_theme = "furo"
html_title = f"bse-recon {release}"
html_theme_options = {
    "source_repository": "https://github.com/Meet-Upadhyay-24/electron-photometric-studio",
    "source_branch": "main",
    "source_directory": "docs/",
}
