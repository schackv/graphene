# -*- coding: utf-8 -*-
"""
Methods and classes specific to representing and calculating lattice properties
from an image.

Created on Tue Sep 24 09:47:05 2013

@author: jsve
"""
import numpy as np
from scipy import ndimage, fftpack
from scipy.spatial.distance import cdist
from skimage.feature import peak_local_max
import cmath

from . import imtools

class parameters:
    theta0 = []
    """Rotation"""
    fourier_extrema = []
    """xy locations of peaks in Fourier power spectrum"""
    fourier_distance = []
    """Avg. distance to peaks in Fourier space"""
    t = []      
    """C-C bond length"""
            
    def compute(self,img,opts):
        """Compute lattice parameters for an image"""

        h,w = img.shape
        if h!=w:
            new_size = np.min((h,w))
            img = img[:new_size,:new_size]

        # Log-power analysis        
        self.PS = powerspectrum(img)
        self.fourier_extrema = _six_local_maxima(self.PS,opts)
              
        # Get radius and rotation 
        radius, angle = _points2avg(self.fourier_extrema, np.array([w*0.5,h*0.5]))
        
        self.fourier_distance = radius
        self.theta0 = angle      
        self.t      = 2*w/(3*radius)
        
      

def powerspectrum(img):
    "Power spectrum of image"
    F = fftpack.fft2(img)
    F = fftpack.fftshift(F)
    ps = np.abs(F)**2   
    return ps      
    
        
def _six_local_maxima(PS,opts):
    """Detect six local maxima in the power spectrum"""
    
    X = ndimage.filters.gaussian_filter(np.log(PS),1)   # Gaussian smoothing
    extrema = peak_local_max(X,**opts)
    
    adjusted_extrema = []
    for i in range(1,7):
        xymin, coef = imtools.fitquadratic(X,extrema[i,],5)
        adjusted_extrema.append(xymin)

    adjusted_extrema = np.vstack(adjusted_extrema)
    adjusted_extrema.shape = 6,2
    
    return adjusted_extrema
        
def _points2avg(points, center):
    """Take the six detected points and estimate radius and rotation in Fourier space"""
        
    npoints = points.shape[0]
    r   = np.empty(npoints)
    phi = np.empty(npoints)
    for i in range(npoints):
        z = complex(points[i,0]-center[0],points[i,1]-center[1])  # "Vectors" from center to peaks
        r[i], phi[i] = cmath.polar(z)   # Polar representation
    
    radius = np.mean(r)     # Avg. radius
    phi[phi<0] += 2*np.pi   # Adjust to [0,2*pi]
    phi = np.sort(phi)  # Sort increasing
    
    phi_diff = [phi[i] - 2*np.pi/6 * i for i in range(npoints) ]
    angle = np.mean(phi_diff)   # Rotation
            
    return radius, angle
        

def hexagonal_centers(im,t):
    """Find hexagonal centers as local minima of im, with an expected hexagonal 
    side length of t.
    
    Returns N x 2 numpy array of coordinates.
    """
    # Contrast enhancement
    
    Y, si = imtools.blob_enhancement(im, np.sqrt(3)*t)
    # Local minima
    nodes = imtools.local_minima(Y, 0.5 * t)
    
    y, x = zip(*nodes)
    
    return np.vstack((x,y)).T, Y, si
        
        
        
        
        
        
        