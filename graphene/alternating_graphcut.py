# -*- coding: utf-8 -*-
"""
Implements the alternating graphcut procedure, which alternates between removing
improbable points and improbable triangles. 

Note that a faster implementation could be obtained by using the max-flow/min-cut
algorithms by Boykov & Kolmogorov, but this requires compilation of a C++ library.
The methods used here are pure python, thanks to the networkx library.

Created on Tue Aug 26 19:39:22 2014

@author: jsve
"""

import networkx as nx
import numpy as np
from scipy.spatial import Delaunay
import copy
import logging
from . import imtools, graphtools


def cleanup(xy,im,alpha=3,beta=2):
    """Alternate between removing improbable simplices and improbable points."""
    
    numchange = np.Inf
    it = 0
    
    # Ensure 1-image is between 0 and 1
    im_neg = imtools.minmax_stretch(1-im,0,1)
    
    xy_hat = copy.deepcopy(xy)
    while numchange > 0:
        N = xy_hat.shape[0]
        dt = Delaunay(xy_hat)
        
        simplices = _remove_improbable_triangles(dt.points,dt.simplices, dt.neighbors, alpha)
        xy_hat = _remove_improbable_points(dt.points,graphtools.tri_edges(simplices),im_neg, beta)
        
        numchange = N - xy_hat.shape[0]
        it+=1
        logging.debug('Alternating graphcut, iteration {}. Removed {} points.'.format(it,numchange))
        
    # Setup graph with these points and simplices
    return xy_hat, simplices
        
        
def _create_graph(source_weights, terminal_weights, nbhood, beta,is_edges=False):
    """Creates a graph with one node per N source_weight, with the given 
    source_weight as capacity from source <--> node. 
    Node <--> node capacities are set to beta for the neighborhood defined by
    the N-long nbhood enumerable."""
    
    D = nx.DiGraph()
    # Set up edges with capacities
    sid = -1        # Id of source
    tid = -2        # Id of sink/terminal
    labels = {}
    for id, (sw, tw, nbs) in enumerate(zip(source_weights, terminal_weights, nbhood)):
        D.add_edge(sid, id, capacity=sw )     # Source weights
        D.add_edge(id, tid, capacity=tw )     # Terminal weights
        labels[id] = '{}'.format(id)
        if not is_edges:
            for nb in nbs:
                if nb<0:  # -1 codes for 'missing neighbor'
                    continue
                D.add_edge(id, nb, capacity=beta)
                
    if is_edges:    # if the neighborhood is specified as edges instead of a list of neighborhood
        for e in nbhood:
            D.add_edge(e[0],e[1], capacity = beta)
            D.add_edge(e[1],e[0], capacity =  beta)
            
    labels[sid] = r'S'
    labels[tid] = r'T'
    
    
    return D
    
        
def _remove_improbable_triangles(xy,tri,nbhood,alpha):
    
    mu = 60
    sigma = 0.85
    
    # Get minimum angle of triangles
    angles = np.rad2deg([minimum_angle(xy[T,:]) for T in tri])

    # Setup Terminal weights and edge weights
    source, terminal = _terminalweights(angles,mu,sigma)
    
    # Cut graph
    D = _create_graph(source,terminal, nbhood, alpha)   
    cut_value, partition = nx.minimum_cut(D,-1,-2)
    reachable, non_reachable = partition
    reachable.discard(-1)   # Remove source
    
    # Keep only the nodes connected to the source    
    tri = tri[list(reachable),:]
    
    return tri
    
def _remove_improbable_points(xy,edges,im, beta):
    
    mu = 1
    sigma=0.15

    # Return intensities under the points
    val = imtools.image_interp(im,xy)
       
    # Setup terminal and edge weights
    source, terminal = _terminalweights(val,mu,sigma)
    
    # Setup graph 
    D = _create_graph(source,terminal, edges, beta, is_edges=True)
    
    # Remove points with no neighbors (i.e., only connected to source and sink)
    deg = D.degree()
    to_remove = [p for p in deg if deg[p]==2 and p>=0]
    D.remove_nodes_from(to_remove)
    
    cut_value, partition = nx.minimum_cut(D,-1,-2)
    reachable, non_reachable = partition
    reachable.discard(-1)   # Remove source
    
    # Return resulting points
    xy = xy[list(reachable),:]
    
    return xy
    

def _terminalweights(x,mu,sigma):
    source_weights = (1-sigma) * np.sqrt(x**2)
    terminal_weights = sigma*np.sqrt((x-mu)**2)
    return source_weights, terminal_weights
    
def _edgeweights(edges,beta):
    return beta
    
    
def minimum_angle(vertices):
    """Calculate minimum angle of a triangle given by 3 vertices"""
    
    min_angle=np.inf
    order = ((0,1,2),(1,0,2),(2,0,1))
    for i in range(3):
        vec1 = vertices[order[i][1],:]-vertices[i,:]
        vec2 = vertices[order[i][2],:]-vertices[i,:]
        angle = angle_between(vec1,vec2)
        min_angle = min(min_angle,angle)
    
    return min_angle

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    angle = np.arccos(np.dot(v1, v2))
    if np.isnan(angle):
        if (v1 == v2).all():
            return 0.0
        else:
            return np.pi
    return angle
    

