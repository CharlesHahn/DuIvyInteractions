# Sphinx 配置 — DuIvyInteractions 中文文档

project = "DuIvyInteractions"
copyright = "2026, DuIvy Team"
author = "DuIvy Team"
version = "0.1.0"
release = "0.1.0"

# -- 扩展 --------------------------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

# -- MyST --------------------------------------------------------------------

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# -- 主题 --------------------------------------------------------------------

html_theme = "sphinx_rtd_theme"

# -- Autodoc -----------------------------------------------------------------

autodoc_member_order = "bysource"
autosummary_generate = True

# -- Intersphinx -------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

# -- 语言 --------------------------------------------------------------------

language = "zh_CN"
