# -- Path setup --------------------------------------------------------------
import os
import sys
sys.path.insert(0, os.path.abspath('../../src/SiFoN'))

# -- Project information -----------------------------------------------------

project = 'SiFoN'
copyright = '2022, Briana Macedo'
author = 'Briana Macedo'

# The full version, including alpha/beta/rc tags
release = '0.0.1'

# -- General configuration ---------------------------------------------------
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon']
templates_path = ['_templates']
exclude_patterns = []
html_theme = 'classic'
html_static_path = ['_static']

napoleon_google_docstring = False
napoleon_numpy_docstring = True
autodoc_dumb_docstring = True
