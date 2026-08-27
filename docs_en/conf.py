# Sphinx configuration — DuIvyInteractions English documentation

project = "DuIvyInteractions"
copyright = "2026, DuIvy Team"
author = "DuIvy Team"
version = "0.1.0"
release = "0.1.0"

# -- Extensions ---------------------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

# -- MyST ---------------------------------------------------------------------

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# -- Theme --------------------------------------------------------------------

html_theme = "sphinx_rtd_theme"

# -- Autodoc ------------------------------------------------------------------

autodoc_member_order = "bysource"
autosummary_generate = True

# -- Intersphinx --------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

# -- Language -----------------------------------------------------------------

language = "en"
