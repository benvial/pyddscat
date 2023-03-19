#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: The Phokaia Developers
# This file is part of pyddscat
# License: GPLv3
# See the documentation at phokaia.gitlab.io/doc/pyddscat

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import pyvista

# necessary when building the sphinx gallery
pyvista.BUILDING_GALLERY = True
pyvista.OFF_SCREEN = True

# Optional - set parameters like theme or window size
pyvista.set_plot_theme("document")

project = "pyddscat"
copyright = "2023, The Phokaia Developers"
author = "The Phokaia Developers"
release = "0.0.1"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.todo",
    "sphinx.ext.graphviz",
    "sphinx.ext.inheritance_diagram",
    "matplotlib.sphinxext.plot_directive",
    "sphinx.ext.napoleon",
    "sphinx_copybutton",
    "sphinx_gallery.gen_gallery",
]


templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_static_path = ["_static"]
html_theme_options = {
    "light_logo": "logo.svg",
    "dark_logo": "logo.svg",
    "sidebar_hide_name": False,
    "light_css_variables": {
        "color-brand-primary": "#df4e38",
        "color-brand-content": "#df4e38",
        "color-sidebar-background": "#ededed",
        "color-sidebar-search-background": "#ededed",
        "font-stack": "IBM Plex Sans Condensed, sans-serif",
        "font-stack--monospace": "IBM Plex Mono, monospace",
    },
    "dark_css_variables": {
        "color-brand-primary": "#df4e38",
        "color-brand-content": "#df4e38",
        "color-sidebar-background": "#2a2d2b",
        "color-sidebar-search-background": "#2a2d2b",
        "color-background-primary": "#3a3b3b",
        "font-stack": "IBM Plex Sans Condensed, sans-serif",
        "font-stack--monospace": "IBM Plex Mono, monospace",
    },
}
html_title = "pyddscat"
html_show_sphinx = False

html_css_files = [
    "css/custom.css",
]
pygments_style = "tango"
pygments_dark_style = "monokai"

sphinx_gallery_conf = {
    "show_signature": False,
    "nested_sections": True,
    # path to your examples scripts
    "examples_dirs": ["../examples"],
    # path where to save gallery generated examples
    "gallery_dirs": ["examples"],
    # directory where function granular galleries are stored
    "backreferences_dir": "generated/backreferences",
    "remove_config_comments": True,
    "filename_pattern": "/plot_",
    "ignore_pattern": r"^((?!/plot_).)*$",  # ignore files that do not start with plot_
    # "ignore_pattern": r"__init__\.py",
    "reference_url": {
        "sphinx_gallery": None,
    },
    "reset_modules": (),
    # "first_notebook_cell": (
    #     "import matplotlib\n" "mpl.style.use('phokaia')\n" "%matplotlib inline"
    # ),
    # "image_scrapers": ("matplotlib", PNGScraper()),
    "image_scrapers": ("pyvista", "matplotlib"),
    # Modules for which function level galleries are created.
    "doc_module": project,
    "thumbnail_size": (800, 800),
    "default_thumb_file": "./_static/logo.png",
    "show_memory": True,
    # "binder": {
    #     "org": "phokaia",
    #     "repo": "phokaia.gitlab.io",
    #     "branch": "master",
    #     "binderhub_url": "https://mybinder.org",
    #     "dependencies": "../environment.yml",
    #     "notebooks_dir": "notebooks",
    #     "use_jupyter_lab": True,
    # },
}
