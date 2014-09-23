# -*- coding: utf-8 -*-
"""
Default settings for various parameters used for gridmatching.

Created on Tue May 28 21:51:43 2013

@author: jsve
"""



_lattice = {'min_distance': 4,   # Size of window for peak detection in FFT space
           'num_peaks': 7}      # Number of peaks to find = 6 + 1 (DC coefficient)

_fine_adjustment = {'beta': 0.7}
          
_annealing = {'T0': 4,
             'Tend': 0.5,
             'reps': 4,
             'nIt': 50}
             
defaults = {'image_border': 20,
            'lattice': _lattice,
            'fine_adjustment': _fine_adjustment,
            'annealing': _annealing}