#!/usr/bin/env python

from distutils.core import setup

setup(
    name='skm-tools',
    version='0.1',
    description='Utilities for leveraging SKM resources',
    author='Carissa Bleker',
    author_email='carissa.bleker@nib.si',
    url='https://github.com/NIB-SI/skm-tools',
    license='GPL3',

    packages=['skm_tools', 'skm_tools.resources'],
    install_requires=[
        'networkX',
        'pandas'
    ],
    extras_require={
        "cytoscape": ["py4cytoscape"],
        "pdf": ["pypdf", "pdfCropMargins", "reportlab"],
    },
)