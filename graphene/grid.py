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
from . import graphtools, misc
import networkx as nx
import logging

class Grid:
    """Defines a Grid base class, consisting of a set of points and a set
    of functions to manipulate these points.
    
    The neighborhood data structure is using the networkx.Graph class."""
    
    def __init__(self, xy,edges=[]):
        self.graph = nx.Graph()
        self.xy = xy
        for id, pos in enumerate(xy):
            self.graph.add_node(id,xy=pos)
        self.graph.add_edges_from(edges)      
       
    @classmethod
    def from_textfile(cls, filename):
        """Grid constructor based on files."""
        xy = np.loadtxt(filename + '.points')
        edgelist = nx.read_edgelist(filename + '.edgelist',nodetype=int)
        return cls(xy,edgelist.edges())

    def write(self, filename):
        """Write the grid to .point and .edgelist files."""
        nx.write_edgelist(self.graph,filename + '.edgelist', data=False)
        np.savetxt(filename + '.points',self.xy)
        
    def resolve_edges(self):
        raise NotImplementedError()
            

    def edges(self):
        return self.graph.edges()            
    
    def edge_lengths(self,scale_xy=(1,1)):
        if self.graph.number_of_nodes()==0 or self.graph.number_of_edges()==0:
            raise NoEdgesException()

        edge_lengths = [misc.weucl(self.graph.node[edge[0]]['xy'], self.graph.node[edge[1]]['xy'],scale_xy) for edge in self.graph.edges_iter()]
#        bondlengths_px = [norm(self.atoms[edge[0],:]-self.atoms[edge[1],:]) for edge in self.atom_edges]
        return np.array(edge_lengths)

    def edge_orientations(self):
        thetas = [misc.orientation(self.graph.node[edge[0]]['xy'], self.graph.node[edge[1]]['xy']) for edge in self.graph.edges_iter()]
        return thetas
    
    """ Add zero-mean Gaussian random noise to the grid points """
    def add_noise(self,noise_std):
        self.xy += np.random.randn(self.xy.shape[0],2)*noise_std
#        for n, attr in self.graph.nodes_iter(data=True):
#            attr['xy'] += np.random.randn(2)*noise_std

    def rotate(self,theta):
        self.xy = rotate_grid(self.xy,theta)
        
    def translate(self,deltaxy):
        self.xy += deltaxy
#        for n, attr in self.graph.nodes_iter(data=True):
#            attr['xy'] += deltaxy

    def plot(self,color='b',linecolor='k',markersize=3):
        import matplotlib.pyplot as plt
        from matplotlib import collections  as mc
        plt.plot(self.xy[:,0],self.xy[:,1],'.',color=color,ms=markersize)
        plt.gca().add_collection(mc.LineCollection(self.line_collection(),colors=linecolor))
        plt.axis('image')

    def line_collection(self):
        return line_collection(self.xy,self.graph.edges())
        
        

class TriangularGrid(Grid):
    """ Implements a triangular grid structure.
    """

    @classmethod
    def from_simplices(cls, xy, simplices):
        edges = graphtools.tri_edges(simplices)
        g1 = cls(xy,edges)
        g1.simplices = simplices
        return g1

    def delaunay_triangulate(self):
        """ Resolve the edges in the current grid leveraging the Delaunay 
        triangulation. Three values of method can be chosen:
        """
        dt = Delaunay(self.xy)
        super().__init__(dt.points,graphtools.tri_edges(dt.simplices))     # init as new graph
        self.simplices = dt.simplices
        return self.edges()
        

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
        super().delaunay_triangulate()
        
        if remove_long:
            # TODO Also remove simplices!
            # Remove too long edges
            lengths = self.edge_lengths()
            idx = lengths > 1.1*np.sqrt(3)*self.t       # Edges to remove
            
            
            E = self.graph.edges()
#            [logging.debug(E[i]) for i in np.where(idx)[0]]
            for i in np.where(idx)[0]:
                self.graph.remove_edge(*E[i])    # Remove edge
            
            
            
            
class HexagonalGrid(Grid):
    
    
    def _simplex_to_arc_nbhood(simplices):
        """Get the simplex-to-arc adjacency matrix as a sparse matrix.
        Input is an iterable of simplices.
        
        This is not pretty, but fairly rapid."""
        
        edges = graphtools.tri_edges(simplices)

        tri_edges = np.sort(edges,axis=1)   # Node with lowest id is first
        simplex_ids = []    # Simplex ids
        edge_ids = []       # Edge ids for each simplex
        for sid, s in enumerate(simplices):
            s_edges = np.sort([[s[0], s[1]], [s[0], s[2]], [s[1], s[2]]],axis=1)    # Edges in triangle
            aux = [np.nonzero((tri_edges[:,0]==edge[0]) & (tri_edges[:,1]==edge[1]))[0] for edge in s_edges]    # Find edge-ids in this particular simplex
            
            edge_ids.append(np.hstack(aux))     # Put in a list
            simplex_ids.append(np.ones(len(aux))*sid)
        simplex_ids = np.hstack(simplex_ids)
        edge_ids = np.hstack(edge_ids)
        
        # Return as sparse matrix
        Ns = len(simplices)
        Ne = len(edges)
        simplex_to_arc = sparse.csc_matrix((np.ones(len(simplex_ids)),np.vstack((simplex_ids,edge_ids))),shape=(Ns,Ne))
        return simplex_to_arc
    
    @classmethod
    def from_triangular(cls,triGrid):
        """Create a hexagonal grid based on its dual, triangular, grid.
        Points are positioned as centers of each simplex.
        """
        
        
        Ne = triGrid.graph.number_of_edges()
        if Ne==0:
            raise NoEdgesException()
            
        logging.debug('Constructing hexagonal grid from triangular ({} edges).'.format(Ne))

        # Get simplex-to-arc neighborhood as sparse matrix
        simplex_to_arc = HexagonalGrid._simplex_to_arc_nbhood(triGrid.simplices)
        
        # Use simplex-to-arc neighborhood to get simplex-to-simplex
        simplex_to_simplex = []
        for eid in range(Ne):
            aux = simplex_to_arc[:,eid].nonzero()[0]
            # Get length-2 combinations (i.e. pairs of simplices)
            [simplex_to_simplex.append(comb) for comb in itertools.combinations(aux, 2)]
        simplex_to_simplex = np.array(simplex_to_simplex)

        # Get simplex-centers        
        xy = []
        for s in triGrid.simplices:
            center = np.mean(triGrid.xy[s,:],axis=0)
            xy.append(center)
        xy = np.vstack(xy)

        obj = cls(xy,edges=simplex_to_simplex)
#        super().__init__(cls,xy=xy,edges=simplex_to_simplex)
        
        logging.debug('Hexagonal grid constructed with {} points'.format(obj.graph.number_of_nodes()))
        return obj
        

        
        
    
def line_collection(xy,edges):
    return np.dstack((xy[list(edges),0],xy[list(edges),1]))

        
        
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
    
    
class NoEdgesException(Exception):
    pass
        

    


    
    
    
    
    
    