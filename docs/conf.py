#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import os
import sys
sys.path.insert(0, os.path.abspath('./..'))
sys.path.append('..')


extensions = ['sphinx.ext.autodoc', 
              'sphinx.ext.napoleon', 
              'sphinx.ext.mathjax', 
              'nbsphinx']
templates_path = ['_templates']

source_suffix = '.rst'
master_doc = 'index'

# General information about the project.
project = 'CrysFiPy'
copyright = '2018, Petr Cermak'
author = 'Petr Cermak'

version = release = '0.5'

language = None

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'CrysFiPydoc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

latex_documents = [
    (master_doc, 'CrysFiPy.tex', 'CrysFiPy Documentation',
     'Petr Cermak', 'manual'),
]


# -- Options for manual page output ---------------------------------------
man_pages = [
    (master_doc, 'crysfipy', 'CrysFiPy Documentation',
     [author], 1)
]


# -- Options for Texinfo output -------------------------------------------
texinfo_documents = [
    (master_doc, 'CrysFiPy', 'CrysFiPy Documentation',
     author, 'CrysFiPy', 'One line description of project.',
     'Miscellaneous'),
]



