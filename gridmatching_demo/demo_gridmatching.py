# -*- coding: utf-8 -*-
"""
This demo shows the entire pipeline of 
1) estimating parameters from the image,
2) finding local minima in the image
3) estimating an initial grid
4) fine adjusting the grid

Created on Fri Jul 25 18:49:53 2014

@author: schackv
"""

import numpy as np
from gridmatching import *
import matplotlib.pyplot as plt
from matplotlib import collections  as mc


def demo_gridmatching():
        
    ## Read DM3 image
    DM3 = io.dm3image('graphene_regular.dm3')
    im = imtools.crop(DM3.image(),601,1100,601,1100)
        
    ## Estimate parameters of lattice
    lp = lattice.parameters()
    lp.compute(im,options.lattice)
    print('Estimated hexagonal side length in pixels = {:8f}'.format(lp.t))
    
    ## Initial points
    extrema, _, _= lattice.hexagonal_centers(im, lp.t)
    
    ## Initial grid
    G = grid.TriangularGrid(extrema)
    G.resolve_edges()    # Delaunay triangulation
    # TODO: Implement actual grid-cleanup-procedure
    print(G.edges)
        
    # Show initial grid
    fig = plt.figure()
    plt.plot(extrema[:,0],extrema[:,1],'.')
    lc = mc.LineCollection(np.dstack((extrema[G.edges,0],extrema[G.edges,1])),colors='k')
    plt.gca().add_collection(lc)
    plt.axis('image')
    plt.title('Initial grid')
    plt.show(block=True)
    
    
    ## Fine adjustment
    cs = bgm.welsh_powell(G.edges)
    # TODO: Color-coding algorithm
    # TODO: 
    
    
    

if __name__=='__main__':
    demo_gridmatching()