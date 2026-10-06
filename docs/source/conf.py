# Configuration file for the Sphinx documentation builder.
#
# For the full list of options see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os

DOCS_DIRECTORY = os.path.dirname(os.path.abspath(__file__))

# -- Project information -----------------------------------------------------

project = 'skm-tools'
copyright = '2023-2026 National Institute of Biology, Slovenia'
author = 'Carissa Bleker'

# The full version, including alpha/beta/rc tags
release = '0.1.0'


# -- General configuration ---------------------------------------------------

extensions = [
    'autoapi.extension',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx_rtd_theme',
    'sphinx.ext.todo',
    'sphinx.ext.intersphinx',

    # embed jupyter notebooks in the documentation
    'nbsphinx',
    'nbsphinx_link',
]

# Autoapi settings
# AutoAPI parses the source without importing it, so the optional dependencies
# (py4cytoscape, pypdf, ...) don't need to be installed to build the docs.
autoapi_dirs = ['../../skm_tools']
# skm_tools/__init__.py imports its submodules; without leaving out 'imported-members',
# AutoAPI would document them a second time under the package.
autoapi_options = [
    'members', 'undoc-members', 'show-inheritance',
    'show-module-summary', 'special-members',
]

# Napoleon settings
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# Intersphinx settings
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'networkx': ('https://networkx.org/documentation/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'py4cytoscape': ('https://py4cytoscape.readthedocs.io/en/latest/', None),
}

# Notebooks are rendered with their saved outputs (many need a running Cytoscape)
nbsphinx_execute = 'never'

# ----------------------------------------------------------------------------

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'style_nav_header_background': '#009739',
}

html_logo = '../figures/logo.png'
html_favicon = '../figures/icon.ico'

html_static_path = ['_static']
html_css_files = ['custom.css']


def ensure_pandoc_installed(_):
    # nbsphinx needs pandoc
    try:
        import pypandoc
        pandoc_dir = os.path.join(DOCS_DIRECTORY, "bin")
        if pandoc_dir not in os.environ["PATH"].split(os.pathsep):
            os.environ["PATH"] += os.pathsep + pandoc_dir
        pypandoc.ensure_pandoc_installed(
            targetfolder=pandoc_dir,
            delete_installer=True,
        )
    except Exception:
        pass  # pandoc already on PATH (CI) or download failed (local SSL) — continue


def setup(app):
    app.connect("builder-inited", ensure_pandoc_installed)
