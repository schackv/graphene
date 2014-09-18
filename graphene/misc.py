# -*- coding: utf-8 -*-
"""
Miscellaneous tools/functions

Created on Thu Sep 18 19:41:13 2014

@author: jsve
"""

import numpy as np
    
def ecdf(x):
    sorted_x=np.sort( x )
    yvals=np.arange(len(sorted_x))/float(len(sorted_x))
    
    return sorted_x, yvals