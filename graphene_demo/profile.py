# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 09:56:23 2014

@author: jsve
"""

import cProfile
import pstats
import demo_simplegrid as dlm


cProfile.run('dlm.demo_simplegrid()','stats.txt')


p = pstats.Stats('stats.txt')
p.strip_dirs().sort_stats(-1).print_stats()

p.sort_stats('cumulative').print_stats(20)