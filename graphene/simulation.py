# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 11:30:28 2014

@author: schackv
"""

import numpy as np
from . import grid
import scipy.ndimage.filters as filters
#PLOT = True

def simulate_image(rows,cols,t,theta=0,sigma_noise=0):
    """ 
    Simulate an image with bright dots in the centers of hexagons with side length t.
    """
    
    # Generate a grid and rotate
    G = grid.SimulatedTriangularGrid(rows,cols,t)
    G.resolve_edges(remove_long=True)
    G.rotate(theta) # Rotate points around center
    G.translate(-np.min(G.xy,axis=0) + t)     # add extra margin
    
    G.add_noise(sigma_noise)
    
    sigma = 0.25*t   # St.d. of gaussian kernel for smoothing the image
    
    # Round for easier image generation
    G.xy = np.round(G.xy-np.min(G.xy,axis=0))   # Ensure that the grid is shifted to positive
    pt_idx = G.xy[:,::-1].astype(int)
    
    # Generate image from centers
    max_xy = np.round(np.max(G.xy,axis=0) * 1.05)   # add five percent
    im = np.zeros((int(max_xy[1]),int(max_xy[0])))    
    im[pt_idx[:,0],pt_idx[:,1]] = 1.0
    
    # Convolve with Gaussian
    im = filters.gaussian_filter(im, sigma)
        
#==============================================================================
#     # Alternative approach
#     xv, yv = np.meshgrid(range(int(max_xy[0])),range(int(max_xy[1])))
#     
#     im = np.zeros(xv.shape)
#     for ctr in G.xy:
#         vals = np.exp(- ( (xv-ctr[0])**2 + (yv-ctr[1])**2) /sigma**2)
#         im += vals
#==============================================================================
    
    
    im = 1-im
    
    return im, G
        

if __name__ == '__main__':
    simulate_image(5,5,10,theta=np.pi*0.25)
    