# -*- coding: utf-8 -*-
"""
This demo shows how to extract global lattice properties 
using the spectral signature in the image.

Created on Fri Jul 25 18:49:53 2014

@author: schackv
"""

import numpy as np
from gridmatching import *
import matplotlib.pyplot as plt




def demo_latticeinfo():
    ## Read DM3 image
    DM3 = io.dm3image('graphene_regular.dm3')
    im = imtools.crop(DM3.image(), 701, 1700, 601, 1600 )   # Use small region to illustrate
    
    # Show image
    fig = plt.figure(figsize=(10,10))    
    plt.matshow(im,fignum=fig.number,cmap=plt.cm.gray)
    plt.show(block=False)
    
    ## Estimate parameters of lattice
    lp = lattice.parameters()
    lp.compute(im,options.lattice)
    print('Hexagonal side length = {:8f}'.format(lp.t))
    print('Rotation in Fourier space = {:4f} radians = {:2f} degrees'.format(lp.theta0, np.rad2deg(lp.theta0)))
        
    # Show lattice power spectrum with found maxima and visualize angle
    h, w = lp.PS.shape
    plt.matshow(np.log(lp.PS),cmap=plt.cm.gray) # Log of power spectrum
    x,y = zip(*lp.fourier_extrema)
    xr = lp.fourier_distance * np.cos(lp.theta0)    
    yr = lp.fourier_distance * np.sin(lp.theta0)
    plt.plot(x,y,'o')
    plt.plot(np.array([0,xr])+w*0.5,np.array([0,yr])+h*0.5,'-xr')
    plt.axis('image')
    plt.show()
    
    

if __name__=='__main__':
    demo_latticeinfo()