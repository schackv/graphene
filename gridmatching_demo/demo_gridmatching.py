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
from scipy import misc
import pickle
import os


def demo_gridmatching():
    
    
    ## Read DM3 image
#    DM3 = io.dm3image('graphene_regular.dm3')
#    im = imtools.crop(DM3.image(),601,1100,601,1100)
    im = misc.imread(r'E:\dtu\phd\graphene\simulated_images\Perfect_dE03eV_Cs-0002mm_df-3mm.png')
    im = imtools.rgb_to_gray(im.astype(np.float32)/255)
    newsize = np.min( (im.shape[0],im.shape[1]) )
    im = im[0:newsize,0:newsize]
    
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
    plt.show(block=False)
    
    
    ## Fine adjustment
    xy_old = G.xy.copy()
    cs = bgm.welsh_powell(G.edges)
    model = bgm.AdaptiveGrid(im, G.edges, G.simplices) 
    xy_hat, E, history = bgm.fit_grid(im,G.xy, G.edges, coding_scheme=cs, beta=0.7, gridenergydefinition=model, opts=options.annealing)
    G.xy = xy_hat
    # Save result to file
    pickle.dump((im,G,history),open('result.pic','wb'))
    
    fig, ax = plt.subplots(1,2)    
    plt.sca(ax[0])
    plt.imshow(im,cmap=plt.cm.gray,interpolation='nearest')
    lc = mc.LineCollection(grid.line_collection(xy_hat,G.edges),colors='k')
    plt.gca().add_collection(lc)
    plt.scatter(xy_old[:,0],xy_old[:,1],c='b',s=30)
    plt.scatter(xy_hat[:,0],xy_hat[:,1],c='r',s=50)
    plt.axis('equal')
    plt.axis('image')
    plt.title('Adjusted grid')
    
    # Show energy over iterations
    plt.sca(ax[1])
    plt.plot(history['energy'],'-')
    plt.xlabel('Iteration')
    plt.ylabel('Energy')
    plt.show(block=False)

    # Place atoms
    G.place_atoms()
    lengths = G.bondlengths_px()
    lengths_nm = lengths * 0.0115 # [nm/pixel]
    print('Bond lengths [px] = {:.6f} +/- {:.6f}'.format(np.mean(lengths),np.std(lengths)))
    
    plt.figure()
    lc = mc.LineCollection(grid.line_collection(G.atoms,G.atom_edges),colors='k')
    plt.gca().add_collection(lc)
    plt.scatter(G.atoms[:,0],G.atoms[:,1],c='b',s=30)
    plt.axis('equal')
    plt.show(block=False)
    
    fix, ax = plt.subplots(1,2)
    plt.sca(ax[0])
    plt.hist(lengths,64)
    plt.sca(ax[1])
    xcdf, F = ecdf(lengths_nm)
    plt.plot(xcdf,F)
    plt.show()
    

    
def ecdf(x):
    sorted_x=np.sort( x )
    yvals=np.arange(len(sorted_x))/float(len(sorted_x))
    
    return sorted_x, yvals

if __name__=='__main__':
    demo_gridmatching()
    place_atoms()

    
    
    