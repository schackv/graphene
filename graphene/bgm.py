# -*- coding: utf-8 -*-
"""
Bayesian grid matching (bgm) module.

Contains methods for minimizing energy in a model of a grid structure.

Created on Sun Jul 27 22:26:56 2014

@author: schackv
"""

import logging
import numpy as np
from copy import deepcopy


class GridEnergyDefinition():
    """Defines a base class for implementations of energy minimization in grids.
    """
    
    def __init__(self, edges, simplices):
        """Initialize the class with a set of edges. Arc-to-arc and node-to-arc 
        neighborhoods are inferred from these edges and stored for later use.
        
        Arc-to-arc neighborhood is defined such that arcs part of the same simplex 
        are neighbors.
        Arc-to-node neighborhood is defined such that edges, where the node is at
        one end, are neighbors to the node.
        """
        edges = np.sort(edges,axis=1)
        self.edges = edges
        
        N = np.max(edges[:])+1
        Narcs = edges.shape[0]

        # Resolve node-to-arc neighborhood
        self.node_to_arc = [np.nonzero(np.any(edges==i,axis=1))[0] for i in range(N)]
        
        # Resolve arc-to-arc neighborhood
        self.arc_to_arc = [set() for i in range(Narcs)]
        for s in simplices:
            s_edges = np.sort([[s[0], s[1]], [s[0], s[2]], [s[1], s[2]]],axis=1)
            aux = [np.nonzero((edges[:,0]==edge[0]) & (edges[:,1]==edge[1]))[0] for edge in s_edges]
            edgeids_in_simplex = np.hstack(aux)     # Edges in this simplex

            # Add to arc sets
            for eid in edgeids_in_simplex:
                self.arc_to_arc[eid] |= set(edgeids_in_simplex)
                
        self.arc_to_arc = [np.array(list(aa)) for aa in self.arc_to_arc]    # Convert to lists
                   
        
        
    def arc_prior(self):
        raise NotImplementedError()
    
    def node_prior(self):
        raise NotImplementedError()
    
    def observation_model(self):
        raise NotImplementedError()
       
       
class AdaptiveGrid(GridEnergyDefinition):
    """Implements an adaptive grid, i.e., a grid that can accommodate 
    spatially changing averages in node-to-node distance.
    """
    
    def __init__(self,im,edges,simplices):
        super().__init__(edges, simplices)
        
        self.im = im
    
    def arc_prior(self,xy):
        """Calculates the mean arc lengths given the node positions xy,
        the arcs and arc_neighborhood defined.
        
        xy                  N x 2
        arcs                Narcs x 2 (indices into xy)
        arc_neighborhood    List of length Narcs with lists of neighbors (indices into arcs).
        """
        # Arc lengths
        arc_lengths = np.sqrt( np.sum( (xy[self.edges[:,0],:] - xy[self.edges[:,1],:])**2,axis=1))
        
        mean_lengths = [np.mean(arc_lengths[nb]) for nb in self.arc_to_arc]

#            
        return np.array(mean_lengths)
        
    def node_prior(self,xy, arc_means, idx=None):
        """Calculates the deviation of each node's distance to its
        neighbors compared to the mean arc length for each arc.
        
        xy      Numpy array of N x 2
        idx     Indices into xy for which this should be returned
        arcs    Numpy array of Narcs x 2 with indices into xy
        node_arc_neighborhood   List of length N with lists of arc neighbors for each node. 
                                Indices are into arcs.
        arc_means   Numpy array of Narcs x 1 with mean arc lengths (as output from arc_prior)

        Returns numpy array of length sum(idx==True) with energies calculated as 
            U_i = 1/n_i * sum_{i~j} (||g_i - g_j|| - \mu_{ij})^2
        """
        
        N = xy.shape[0]
        if idx is None:     # Not given, return for all nodes
            idx = np.array(range(N))
            
        Nidx = len(idx)
        node_arc_neighborhood = [self.node_to_arc[i] for i in idx]
        
        # Arc lengths
        arc_lengths = np.sqrt( np.sum((xy[self.edges[:,0],:] - xy[self.edges[:,1],:])**2,axis=1))
        
        # Squared difference from mean
        Dsq = (arc_lengths - arc_means)**2
        
        # Squared difference from mean for each neighbor and average
        U = [np.mean(Dsq[nb]) for nb in node_arc_neighborhood]
        return np.array(U)
        
        
    def observation_model(self, xy):
        """Calculates the energy in the image at the given positions."""
        m, n = self.im.shape
        
        xy = np.fmin(xy, [n-1, m-1])
        xy = np.fmax(xy, [0,0])
        
        values = [self.im[y,x] for x,y in np.round(xy)]
        
        return np.array(values)


class FixedGrid(AdaptiveGrid):
    """Implements a fixed grid, i.e., a grid that expects a the same distance
    between nodes throughout the grid.
    
    This class inherits and leverages the behavior of AdaptiveGrid, but rather than
    calculating the expected arc length around each arc, the fixed length is 
    returned.
    """

    def __init__(self,im,expected_distance,edges,simplices):
        super().__init__(im, edges, simplices)
        self.expected_distance = np.array(expected_distance)
        
    def arc_prior(self,xy):
        """Not used in this implementation."""
        return self.expected_distance
    
    

def take_step(xy):
    """Take a random step in one of eight directions for each point in xy."""
    
    N = xy.shape[0]
    steps = ([1,0], [-1,0], [0,-1],[0,1],[1,1],[-1,1],[-1,-1],[1,-1])
    nsteps = len(steps)
    
    rndint = np.random.randint(0,nsteps,N)
    xy_new = np.array([pt+steps[i] for pt, i in zip(xy,rndint)])
    return xy_new

def simulated_annealing(xy, coding_scheme, beta, model, T0=4, Tend=0.5, reps=4,nIt=500,max_nz=5):
    """Simulated annealing of the N sites in xy divided into coding_scheme.
    The regularization parameter beta is used to weigh observation versus geometry.
    
    model is an instance of the GridEnergyDefinition class
    defining the energy functions to be minimized.
    
    Additional parameters include 
    T0      start temperature 
    Tend    final temperature
    reps    number of repetitions on each temperature
    nIt     maximum number of iterations
    max_nz  maximum number of iterations with zero accepted moves, before stopping
    """
    
    # Number of colors in coding scheme 
    coding_scheme = coding_scheme
    ncodes = np.max(coding_scheme)
    
    # Temperature lowering constant
    C = (Tend/T0)**(1/nIt)
    
    # Get initial energies
    pA = model.arc_prior(xy)
    uN = model.node_prior(xy,pA)
    uObs = model.observation_model(xy)
    E = beta*uN + (1-beta)*uObs
    
    ## Optimization
    T = T0
    it = 0
    nz = 0      # Number of zero-moves in a row
    n_accepted = np.inf
    history = {'xy': [], 'energy': [], 'arcprior': []}
    
    while (T>Tend) & (it<nIt) & (nz<max_nz):
        it += 1
        
        n_accepted = 0
        for r in range(reps):   # Repeat on each temperature
            for c in range(ncodes+1):   # For each color group
                idx = np.where(coding_scheme==c)[0]
                ngroup = len(idx)
                
                # Take a new step for each node
                xy_new = np.copy(xy)
                xy_new[idx,:] = take_step(xy_new[idx,:])
                
                # Get geometry and observation energy
                uN = model.node_prior(xy_new,pA,idx=idx)
                uObs = model.observation_model(xy_new[idx,:])
                
                # Energy differences at these nodes
                E_new = beta*uN + (1-beta)*uObs     # this is ngroup x 1
                dE = E_new - E[idx]
                
                # Accept the step?
                p = np.fmin(1,np.exp(-1/T * dE)) # min(1,p)
                accept = p >= np.random.uniform(size=ngroup)
                accept = dE < 0
                n_accepted += np.sum(accept)
                
                # Update configuration
                xy[idx[accept],:] = xy_new[idx[accept],:]
                E[idx[accept]] = E_new[accept]  # Transfer calculated energies
                
            
        # Save to history for later visualization
        history['xy'].append(xy.copy())
        history['energy'].append(np.sum(E))
        history['arcprior'].append(pA.copy())
            
        logging.debug('Iteration {:g}. Number accepted = {:g}, energy = {:.8f}'.format(it, n_accepted,np.mean(E)))
        
        if n_accepted==0: 
            nz+=1
        else:
            nz=0
        
        # Update arc prior
        pA = model.arc_prior(xy)
                
        # Decrease temperature
        T *= C
    
    return xy, E, history
        

    
def _edges_to_neighborlists(edges):
    
    N = np.max(edges[:])+1  # Number of nodes
    neighbor_lists = [[] for i in range(N)]
    for e in edges:
        neighbor_lists[e[0]].append(e[1])
        neighbor_lists[e[1]].append(e[0])
        
    return neighbor_lists
    
def welsh_powell(edges):
    """The Welsh-Powell algorithm for setting up a coding scheme
    for the nodes referred to by edges. 
    
    edges       Narcs x 2
    
    Returns a np.max(edges[:]) vector of integers from 0 to the necessary number
    of colors (minimum is the chromatic number of the graph).
    """
    
    neighbor_lists = _edges_to_neighborlists(edges)
    N = len(neighbor_lists)
    
    degree = np.array([len(nblist) for nblist in neighbor_lists])
    grouping = np.empty(N,dtype=int)
    grouping[:] = np.inf
    
    sort_order = degree.argsort()[::-1] # Sort descending
    for idx in sort_order:
        used_colors = grouping[neighbor_lists[idx]]
        for j in range(N):  # Find first color not used
            if not np.any(used_colors==j):
                grouping[idx] = j
                break
    return grouping
    
    
def fit_grid(im,xy,edges,coding_scheme=None,beta=0.5,gridenergydefinition=None, anneal_opts={}):
    """Fit the grid defined by xy and edges to the image im.
    
    im is an m x n image matrix
    xy is a numpy N x 2 array
    edges is a numpy N_edges x 2 array, indexing into xy
    coding_scheme is an N long iterable with integers indicating colors of nodes. 
        If coding_scheme = None, it is estimated using the Welsh-Powell algorithm
    
    The gridenergydefinition implements the energy functions defining the model. 
    These functions are called during simulated annealing.
    If gridenergydefinition=None the standard model is chosen.
    """
    
    im = (im-np.mean(im))/(np.std(im))
    
    if coding_scheme is None:
        coding_scheme = welsh_powell(edges)
        
    if gridenergydefinition is None:
        gridenergydefinition = AdaptiveGrid(im,edges)
        
    xy_out, E, history = simulated_annealing(deepcopy(xy),coding_scheme,beta=beta,model=gridenergydefinition, **anneal_opts)
    return xy_out, E, history
    
    
    