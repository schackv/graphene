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
from graphene import *
import matplotlib.pyplot as plt

from matplotlib import collections  as mc

import pickle


def demo_gridmatching():
    
    filename = 'graphene_regular.dm3'
#    filename = r'E:\dtu\phd\graphene\simulated_images\Perfect-2_dE03eV_Cs-0002mm_df-3mm_0142.png'
    ## Read DM3 image
    im = imtools.read_image(filename)
    im = imtools.crop(im,601,1100,601,1100)
    
    # Make quadratic
    newsize = np.min( (im.shape[0],im.shape[1]) )
    im = im[0:newsize,0:newsize]
    
    ## Estimate parameters of lattice
    lp = lattice.parameters()
    lp.compute(im,options.lattice)
    print('Estimated hexagonal side length in pixels = {:8f}'.format(lp.t))
    
    ## Initial points
    extrema, _, _= lattice.hexagonal_centers(im, lp.t)
    
    ## Initial grid
    xy, simplices = alternating_graphcut.cleanup(extrema,im)
    G = grid.TriangularGrid.from_simplices(xy,simplices)
        
    # Show initial grid
    fig = plt.figure()
    plt.plot(extrema[:,0],extrema[:,1],'.k',label='Extrema')
    plt.plot(G.xy[:,0],G.xy[:,1],'ob',label='Maintained points')
    plt.gca().add_collection(mc.LineCollection(G.line_collection(),colors='k'))
    plt.axis('image')
    plt.title('Initial grid')
    plt.legend()
    plt.show(block=False)
    
    
    ## Fine adjustment
    xy_old = G.xy.copy()
    cs = bgm.welsh_powell(G.edges())
    model = bgm.AdaptiveGrid(im, G.edges(), G.simplices) 
    xy_hat, E, history = bgm.fit_grid(im,G.xy, G.edges(), coding_scheme=cs, beta=0.7, gridenergydefinition=model, opts=options.annealing)
    G.xy = xy_hat
    # Save result to file
    pickle.dump((im,G,history),open('result.pic','wb'))
    
    fig, ax = plt.subplots(1,2)    
    plt.sca(ax[0])
    plt.imshow(im,cmap=plt.cm.gray,interpolation='nearest')
    lc = mc.LineCollection(G.line_collection(),colors='k')
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


#    im, G, history = pickle.load(open('result.pic','rb'))

    # Place atoms
    H = grid.HexagonalGrid.from_triangular(G)
    lengths = H.edge_lengths()
    lengths_nm = lengths * 0.0115 # [nm/pixel]
    print('Bond lengths [px] = {:.6f} +/- {:.6f}'.format(np.mean(lengths),np.std(lengths)))
    
    plt.figure()
    lc = mc.LineCollection(H.line_collection(),colors='k')
    plt.gca().add_collection(lc)
    plt.scatter(H.xy[:,0],H.xy[:,1],c='b',s=30)
    plt.axis('equal')
    plt.show(block=False)
        
    fix, ax = plt.subplots(1,2)
    plt.sca(ax[0])
    plt.hist(lengths,64)
    plt.sca(ax[1])
    xcdf, F = misc.ecdf(lengths_nm)
    plt.plot(xcdf,F)
    plt.show(block=True)
    


if __name__=='__main__':
    demo_gridmatching()


    
    
    