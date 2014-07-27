# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 11:30:28 2014

@author: schackv
"""

import numpy as np
import grid
PLOT = True
""" 
Simulate an image with bright dots in the centers of hexagons with side length t.
"""
def simulate_image(rows,cols,t,theta=0,sigma_noise=0):
    
    # Generate a grid and rotate
    G = grid.TriangularGrid(rows,cols,t)
    G.rotate(theta) # Rotate points around center
    G.translate(-np.min(G.xy,axis=0) + t)     # add extra margin
    G.add_noise(sigma_noise)
    E = G.resolve_edges()
    sigma = 0.25*t   # St.d. of gaussian kernel for smoothing the image
        
    # Generate image from centers
    max_xy = np.round(np.max(G.xy,axis=0) * 1.05)   # add five percent
    xv, yv = np.meshgrid(range(int(max_xy[0])),range(int(max_xy[1])))
    
    im = np.zeros(xv.shape)
    for ctr in G.xy:
        vals = np.exp(- ( (xv-ctr[0])**2 + (yv-ctr[1])**2) /sigma**2)
        im += vals        
        

if __name__ == '__main__':
    simulate_image(5,5,10,theta=np.pi*0.25)
    