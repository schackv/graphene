# -*- coding: utf-8 -*-
"""
Miscellaneous tools/functions

Created on Thu Sep 18 19:41:13 2014

@author: jsve
"""

import numpy as np
import json
    
def ecdf(x):
    sorted_x=np.sort( x )
    yvals=np.arange(len(sorted_x))/float(len(sorted_x))
    
    return sorted_x, yvals
    
def weucl(xy1,xy2,weight_xy=(1,1)):
    """Weighted Euclidean distance between xy1 and xy2."""
    
    return np.sqrt( np.sum( ((xy1-xy2)*weight_xy)**2))
    
def orientation(xy1,xy2):
    """Orientation of the vector from xy1 to xy2."""
    vec = xy2-xy1
    theta = np.arctan2(vec[1],vec[0])
    return theta
    
def angular_diff(theta1, theta2):
    """Smallest angular difference between theta1 and theta2."""
    return np.angle(np.exp(1j*theta1) / np.exp(1j*theta2))
    
    
def circular_binning(angles,nbins=6,theta0=np.pi/6,half_circle=False):
    """Bin the angles in nbins bins, equidistantly spaced around the circle 
    starting at theta0.
    
    If half_circle=True, the returned number of bins will be halved, i.e.,
    only the first nbins/2 bins will be returned. 
    nbins needs to be even in this case.
    
    Returns bin_centers, bin_idx"""
    
    delta = 2*np.pi/nbins
    bin_centers = np.arange(nbins)*delta + theta0
    bin_idx = [np.argmin( np.abs(angular_diff(angle,bin_centers))) for angle in angles]
    
    if half_circle:
        bin_centers = bin_centers[0:nbins/2]
        new_items = [x if x < nbins/2 else x-nbins/2 for x in bin_idx]
        bin_idx = new_items
        
    return bin_centers, bin_idx


def _readdict(filename):
    """Read a dictionary from file using JSON."""
    with open(filename, 'r') as f:
        d = json.load(f)
    return d

def _writedict(filename,config_dict):
    """Write a dictionary to file using JSON."""
    
    with open(filename,'w') as f:
        json.dump(config_dict,f,indent=4)
        
def _merge(a,b,path=None):
    "merges b into a"
    if path is None: path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                _merge(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass # same leaf value
            else:   # a values are overwritten by b values
                a[key] = b[key]
#                raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a