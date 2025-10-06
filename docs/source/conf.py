# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

FILE_HEADER = """Copyright (C) 2024 Jessie Fielding

This file is part of simble.

simble is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

simble is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with simble.  If not, see <https://www.gnu.org/licenses/>."""

import os
import sys
sys.path.insert(0, os.path.abspath('../../../simble/'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'simble'
copyright = '2025, Jessie Fielding'
author = 'Jessie Fielding'
release = '0.0.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        'sphinx.ext.autodoc',
        'sphinx.ext.apidoc',
        'sphinx.ext.napoleon', # Uncomment if using Napoleon for docstring styles
    ]

templates_path = ['_templates']
exclude_patterns = ['/data', '/dev_helper.py']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = [
        'custom.css',
    ]
html_sidebars = { '**': ['globaltoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html'] }


def remove_initial_string(app, what, name, obj, options, lines):
    if what == "module":
        if lines:
            joined_lines = "\n".join(lines)
            if FILE_HEADER.strip() in joined_lines.strip():
                joined_cleaned = joined_lines.replace(FILE_HEADER.strip(), "")
                lines.clear()
                lines.extend(joined_cleaned.strip().split("\n"))
                print("File header found")
    else:
        if lines and "\n".join(lines).strip() == FILE_HEADER.strip():
            print("File header also found")
            return []

def setup(app):
    app.connect("autodoc-process-docstring", remove_initial_string)
