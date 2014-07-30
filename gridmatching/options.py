# -*- coding: utf-8 -*-
"""
Default settings for various parameters used for gridmatching.

Created on Tue May 28 21:51:43 2013

@author: jsve
"""



lattice = {'min_distance': 4,   # Size of window for peak detection in FFT space
           'num_peaks': 7}      # Number of peaks to find = 6 + 1 (DC coefficient)
           
annealing = {'T0': 4,
             'Tend': 0.5,
             'reps': 4,
             'nIt': 100}