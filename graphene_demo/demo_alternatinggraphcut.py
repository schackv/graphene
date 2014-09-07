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



def demo_alternatinggraphcut():
    t = 12
    ## Simulate small image
    im, G = simulation.simulate_image(10,20,t,theta=0.4*np.pi )
    
    # Show grid overlaid image
    plt.imshow(im,cmap=plt.cm.gray)
    plotgrid(G)
    plt.show()
    
    # Add random points
    Nnew = 100
    x_new = np.random.uniform(np.min(G.xy[:,0]),np.max(G.xy[:,0]), Nnew)
    y_new = np.random.uniform(np.min(G.xy[:,1]),np.max(G.xy[:,1]), Nnew )
    xy_new = np.vstack((x_new,y_new)).T
    
    # Add to grid and make naive guess for edges
    G.xy = np.vstack((G.xy,xy_new))
    G.resolve_edges(remove_long=False)
    
    # Show this modified grid
    plotgrid(G)
    plt.show()
    
    # Clean up using alternating graph cut
    xy_hat, simplices_hat = graph.alternating_graphcut(G.xy,im,alpha=3,beta=2)
    G_hat = grid.Grid.from_simplices(xy_hat,simplices_hat)
    plotgrid(G_hat)
    plt.show()
    
#    print(xy_hat)
#   print(simplices_hat)
    


def plotgrid(G):
    plt.plot(G.xy[:,0],G.xy[:,1],'.b')
    lc = mc.LineCollection(G.line_collection(),colors='k')
    plt.gca().add_collection(lc)
   
    

if __name__=='__main__':
    demo_alternatinggraphcut()