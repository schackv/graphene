# -*- coding: utf-8 -*-
"""
This demo shows the steps in finding local minima in the image

Created on Sun Jul 25 18:49:53 2014

@author: schackv
"""

from graphene import *
import matplotlib.pyplot as plt


def demo_localminima():
        
    opts = options.defaults
    
    ## Read DM3 image
    im = imtools.read_image('graphene_regular.dm3')
    im = imtools.crop(im,601,1100,601,1100)
        
    ## Estimate parameters of lattice
    lp = lattice.parameters()
    lp.compute(im,opts['lattice'])
    print('Estimated hexagonal side length in pixels = {:8f}'.format(lp.t))
    
    # Get local minima  
    extrema, contrast_enhanced, si = lattice.hexagonal_centers(im, lp.t)

    print('{:g} initial points found'.format(len(extrema)))
    print(extrema)
    
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
    plt.plot(extrema[:,0],extrema[:,1],'o',ms=3)
    plt.axis('image')
    plt.show(block=True)
    

if __name__=='__main__':
    demo_localminima()