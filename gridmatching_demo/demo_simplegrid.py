# -*- coding: utf-8 -*-
"""
This demo shows 
1) Simulating a small image with responses in a triangular pattern
2) fitting a grid using bayesian grid matching

Created on Fri Jul 25 18:49:53 2014

@author: schackv
"""

import numpy as np
from gridmatching import *
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from scipy.linalg import norm


def demo_simplegrid():
        
   
    t = 12
    ## Simulate small image
    im, G = simulation.simulate_image(10,5,t,theta=0.4*np.pi )
    
    # Show image with true grid overlaid
    plt.figure()
    plt.imshow(im,cmap=plt.cm.gray)
    plt.plot(G.xy[:,0],G.xy[:,1],'.b')
    lc = mc.LineCollection(G.line_collection(),colors='k')
    plt.gca().add_collection(lc)
    plt.colorbar()
    plt.axis('equal')
    plt.axis('image')
    plt.show(block=False)

    # Translate the grid and add noise
    true_xy = G.xy.copy()
    G.translate(np.array([3,4]))
    G.add_noise(1.5)
    
    # Coding scheme
    cs = bgm.welsh_powell(G.edges)
    fig, ax = plt.subplots(2,2)
    plt.sca(ax[0,0])
    lc = mc.LineCollection(G.line_collection(),colors='k')
    plt.gca().add_collection(lc)
    plt.scatter(G.xy[:,0],G.xy[:,1],c=cs,s=150)
    plt.title('Coding scheme')
    plt.axis('equal')
    plt.axis('image')
    
    # Setup an adaptive grid model
    model = bgm.AdaptiveGrid(im, G.edges, G.simplices) 
    xy_hat, E, history = bgm.fit_grid(im,G.xy, G.edges, coding_scheme=cs, beta=0, gridenergydefinition=model, opts=options.annealing )
    
    [print('{:.4f}'.format( norm(x1-x2) ) ) for x1, x2 in zip(xy_hat[G.edges[:,0]],xy_hat[G.edges[:,1]])]
    
    plt.sca(ax[0,1])
    plt.imshow(im,cmap=plt.cm.gray,interpolation='nearest')
    lc = mc.LineCollection(grid.line_collection(xy_hat,G.edges),colors='k')
    plt.gca().add_collection(lc)
    plt.scatter(G.xy[:,0],G.xy[:,1],c='b',s=30)
    plt.scatter(xy_hat[:,0],xy_hat[:,1],c='r',s=50)    
    plt.scatter(true_xy[:,0],true_xy[:,1],c='w',s=30)
    plt.axis('equal')
    plt.axis('image')
    plt.title('Adjusted grid')
    
    # Show energy over iterations
    plt.sca(ax[1,0])
    plt.plot(history['energy'],'-')
    plt.xlabel('Iteration')
    plt.ylabel('Energy')
    
    # Show distance to true xy
    plt.sca(ax[1,1])
    plt.plot([np.mean(np.abs(xy-true_xy)) for xy in history['xy']])
    plt.xlabel('Iteration')
    plt.ylabel('Average absolute distance to true positions')
    plt.show(block=True)
  
    
    
    

if __name__=='__main__':
    demo_simplegrid()