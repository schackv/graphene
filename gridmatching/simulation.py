# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 11:30:28 2014

@author: schackv
"""

import numpy as np
PLOT = True
""" 
Simulate an image with bright dots in the centers of hexagons with side length t.
"""
def simulate_image(rows,cols,t,theta=0):
    sigma = 0.25*t   # St.d. of gaussian kernel for smoothing
    
    # Generate centers
    xy = []
    for i in range(rows):
        for j in range(cols):
            xy.append(position(i,j,t))
    xy = np.vstack(xy)
    
    # Rotate points around center
    xy = rotate_grid(xy,theta)
    xy -= np.min(xy,axis=0) - t     # add extra margin
    if PLOT:
        import matplotlib.pyplot as plt
        plt.plot(xy[:,0],xy[:,1],'x')
        plt.axis('equal')
#        plt.show()
        
    # Generate image from centers
    max_xy = np.round(np.max(xy,axis=0) * 1.05)   # add five percent
    xv, yv = np.meshgrid(range(int(max_xy[0])),range(int(max_xy[1])))
    
    im = np.zeros(xv.shape)
    for ctr in xy:
        vals = np.exp(- ( (xv-ctr[0])**2 + (yv-ctr[1])**2) /sigma**2)
        im += vals        
        
    if PLOT:
        plt.imshow(im,cmap='gray')
        plt.axis('image')
        plt.show()
    
    
""" Rotate the grid around its center """
def rotate_grid(xy,theta):
    ctr = np.mean(xy,axis=0)
    R = np.matrix([[np.cos(theta),-np.sin(theta)],
                   [np.sin(theta),np.cos(theta)]])

    xy_rotated = (xy-ctr) * R.T
    return np.array(xy_rotated + ctr)     # readd center

""" Get the position of the hexagon center at a given row and column idx"""
def position(row_idx,col_idx,t):
    x = t + col_idx*t
    y = tri_sidelength(t) * row_idx
    if col_idx % 2 == 0:    # Even rows
        y += 0.5*tri_sidelength(t)
    else:
        y += tri_sidelength(t)
    
    
    return x, y

"""Side length of triangles connecting centers of hexagons with side length t"""
def tri_sidelength(t):
    return np.sqrt(3)*t
    

if __name__ == '__main__':
    simulate_image(5,5,10,theta=np.pi*0.25)