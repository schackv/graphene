# -*- coding: utf-8 -*-
"""
Misc graph methods

Created on Thu Sep 11 08:50:48 2014

@author: jsve
"""

import numpy as np

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