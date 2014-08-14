# -*- coding: utf-8 -*-
"""
Contains abstractions and implementations of a grid structure.

Created on Wed Jun 25 13:33:34 2014

@author: schackv
"""
import numpy as np
from scipy.spatial import Delaunay
from scipy.linalg import norm
from scipy import sparse
import itertools

"""Defines a Grid base class, consisting of a set of points and a set
of functions to manipulate these points"""
class Grid:
    
    def __init__(self, xy,edges=[]):
        self.xy = xy
        self.edges = edges     # Initialize edges as empty
        self.atoms = []
        self.atom_edges = []
        
    def resolve_edges(self):
        raise NotImplementedError()
        
    def place_atoms(self):
        if len(self.edges)==0:
            raise NoEdgesException()
        
        # Get simplex-to-arc neighborhood as sparse matrix
        self.edges = np.sort(self.edges,axis=1)
        simplex_ids = []
        edge_ids = []
        for sid, s in enumerate(self.simplices):
            s_edges = np.sort([[s[0], s[1]], [s[0], s[2]], [s[1], s[2]]],axis=1)
            aux = [np.nonzero((self.edges[:,0]==edge[0]) & (self.edges[:,1]==edge[1]))[0] for edge in s_edges]
            aux = np.hstack(aux)
            edge_ids.append(aux)
            simplex_ids.append(np.ones(len(aux))*sid)
        simplex_ids = np.hstack(simplex_ids)
        edge_ids = np.hstack(edge_ids)
        
        Ns = len(self.simplices)
        Ne = len(self.edges)
        simplex_to_arc = sparse.csc_matrix((np.ones(len(simplex_ids)),np.vstack((simplex_ids,edge_ids))),shape=(Ns,Ne))

        # Use simplex-to-arc neighborhood to get simplex-to-simplex
        simplex_to_simplex = []
        for eid in range(Ne):
            aux = simplex_to_arc[:,eid].nonzero()[0]
            # Get length-2 combinations (i.e. pairs of simplices)
            [simplex_to_simplex.append(comb) for comb in itertools.combinations(aux, 2)]
        simplex_to_simplex = np.array(simplex_to_simplex)
        
        # One atom per simplex!
        self.atoms = []
        for s in self.simplices:
            atom = np.mean(self.xy[s,:],axis=0)
            self.atoms.append(atom)
        self.atoms = np.vstack(self.atoms)
        self.atom_edges = simplex_to_simplex # Same neighborhood as simplices
            
    def bondlengths_px(self):
        if len(self.atoms)==0 or len(self.atom_edges)==0:
            raise NoEdgesException()
        
        bondlengths_px = [norm(self.atoms[edge[0],:]-self.atoms[edge[1],:]) for edge in self.atom_edges]
        return bondlengths_px

    
    """ Add zero-mean Gaussian random noise to the grid points """
    def add_noise(self,noise_std):
        self.xy += np.random.randn(*self.xy.shape)*noise_std

    def rotate(self,theta):
        self.xy = rotate_grid(self.xy,theta)
        
    def translate(self,deltaxy):
        self.xy += deltaxy

    def line_collection(self):
        return line_collection(self.xy,self.edges)

class NoEdgesException(Exception):
    pass
        
class TriangularGrid(Grid):
    """ Implements a triangular grid structure.
    """
        
    def resolve_edges(self, remove_long=True):
        """ Resolve the edges in the current grid using the Delaunay triangulation.
        """
        dt = Delaunay(self.xy)
        
        self.edges = tri_edges(dt.simplices)
        self.simplices = dt.simplices
        return self.edges
        

class SimulatedTriangularGrid(TriangularGrid):
    """Represents a triangular grid with a given number of rows and columns.
    
    Inherits TriangularGrid.
    """
    def __init__(self, rows,cols, t):
        self.t = t      # Hexagonal side length
        # Generate centers
        xy = []
        for i in range(rows):
            for j in range(cols):
                xy.append(center_position(i,j,t))
        xy = np.vstack(xy)
        super().__init__(xy)
        
    def resolve_edges(self,remove_long=True):
        super().resolve_edges()
        
        if remove_long:
            # Remove too long edges
            E = self.edges
            xy1 = self.xy[E[:,0],:]
            xy2 = self.xy[E[:,1],:]
            lengths = [norm(x1-x2) for x1, x2 in zip(xy1,xy2)]
            idx = lengths < 1.1*np.sqrt(3)*self.t
            self.edges = E[idx,:]
            
def line_collection(xy,edges):
    return np.dstack((xy[list(edges),0],xy[list(edges),1]))

def tri_edges(simplices):
    """Return a list of unique edges representing the simplices."""
    edges = set()

    def add_edge(i, j):
        """Add a line between the i-th and j-th points, if not in the list already"""
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add( (i, j) )

    # loop over triangles: 
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in simplices:
        add_edge(ia, ib)
        add_edge(ib, ic)
        add_edge(ic, ia)
    
    return np.array(list(edges))
        
        
""" Get the position of the hexagon center at a given row and column idx"""
def center_position(row_idx,col_idx,t):
    x = 1.5* t*(1 + col_idx)
    y = tri_sidelength(t) * row_idx
    if col_idx % 2 == 0:    # Even rows
        y += 0.5*tri_sidelength(t)
    else:
        y += tri_sidelength(t)
    
    return x, y
        
"""Side length of triangles connecting centers of hexagons with side length t"""
def tri_sidelength(t):
    return np.sqrt(3)*t
    
        
""" Rotate a set of points around their center """
def rotate_grid(xy,theta):
    ctr = np.mean(xy,axis=0)
    R = np.matrix([[np.cos(theta),-np.sin(theta)],
                   [np.sin(theta),np.cos(theta)]])

    xy_rotated = (xy-ctr) * R.T
    return np.array(xy_rotated + ctr)     # readd center