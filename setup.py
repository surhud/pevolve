#!/usr/bin/env python

"""
setup.py file for pevolve
"""

from distutils.core import setup, Extension

cosmology_module = Extension('_cosmology',
                           sources=['src/cosmology.i', 'src/cosmology.cpp'],
                           swig_opts=["-c++"],
                           libraries=['m','gsl','gslcblas'],
                           )

setup (name        = 'pevolve',
       version     = '1.0rc',
       author      = "Surhud More",
       url         = "http://member.ipmu.jp/surhud.more/research",
       author_email= "surhud.more@ipmu.jp",
       description = """Pseudo-evolution and cosmology module""",
       ext_modules = [cosmology_module],
       license     = ['GPL'],
       py_modules  = ["cosmology"],
       package_dir = { '':'src'},
       )

