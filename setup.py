#!/usr/bin/env python

import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

print(find_packages())
setup(
    name = 'graphene',
    version = '0.1',
    author = 'Jacob Schack Vestergaard',
    author_email = 'jsve@dtu.dk',
    description = "Tools to fit grid models to graphene images.",
    license = 'MIT',
    url = 'http://compute.dtu.dk/~jsve',
    packages = find_packages(),
    install_requires = ['numpy', 'scipy', 'matplotlib', 'dm3lib_v099','skimage','networkx'],
    long_description = read('README.md'),
    classifiers = [
        'Development Status :: 1 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
    ]
)
