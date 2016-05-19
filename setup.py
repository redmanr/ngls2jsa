# -*- coding: utf-8 -*-
"""
Created on Wed May 18 19:12:12 2016

@author: russell
"""

from setuptools import setup, find_packages
import sys

if sys.version_info[0] < 3:
    print('python version: ' + sys.version_info)
    print('The jpggps2kml package is only compatible with Python version 3.n')
    sys.exit(-1)

setup(name='ngls2jsa',
      version='0.1.0',
      description='scripts to decorate NGLS files with FITS headers needed '
                  'for ingestion into the JCMT Science Archive,
      author='Russell O. Redman',
      author_email='russell.o.redman@gmail.com',
      install_requires=['tools4caom2',
                        'jcmt2caom2'],
      packages=find_packages(exclude=['*.test']),
    )