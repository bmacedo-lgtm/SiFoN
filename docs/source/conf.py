import os
import sys
sys.path.insert(0, os.path.abspath('../../src/SiFoN'))

# -- Project information -----------------------------------------------------

project = 'SiFoN'
copyright = '2022, Briana Macedo'
author = 'Briana Macedo'
release = '0.0.1'


# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    'sphinx.ext.napoleon'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
language = 'en'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

html_theme = 'press'
# html_static_path = ['_static']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
todo_include_todos = True