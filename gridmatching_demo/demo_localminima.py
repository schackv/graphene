# -*- coding: utf-8 -*-
"""
This demo shows the steps in finding local minima in the image

Created on Sun Jul 25 18:49:53 2014

@author: schackv
"""

import numpy as np
from gridmatching import *
import matplotlib.pyplot as plt


def demo_localminima():
        
    ## Read DM3 image
    DM3 = io.dm3image('graphene_regular.dm3')
    im = imtools.crop(DM3.image(),601,1100,601,1100)
        
    ## Estimate parameters of lattice
    lp = lattice.parameters()
    lp.compute(im,options.lattice)
    print('Estimated hexagonal side length in pixels = {:8f}'.format(lp.t))
    
    # Get local minima  
    extrema, contrast_enhanced, si = lattice.hexagonal_centers(im, lp.t)

    # Show steps
    fig, ax = plt.subplots(nrows=1,ncols=3)
    plt.sca(ax[0])
    plt.matshow(si,fignum=False, cmap=plt.cm.gray)   
    plt.title('Shape index (blob enhancer)')
    plt.sca(ax[1])
    plt.matshow(contrast_enhanced,fignum=False, cmap=plt.cm.gray)
    plt.title('Contrast enhanced')
    plt.sca(ax[2])
    plt.matshow(im, fignum=False, cmap=plt.cm.gray)
    y,x = zip(*extrema)
    plt.plot(x,y,'o',ms=3)
    plt.axis('image')
    plt.show()
    
    print(extrema)
    
    
    
    

if __name__=='__main__':
    demo_localminima()