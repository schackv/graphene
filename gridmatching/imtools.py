# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 13:43:54 2013

@author: jsve
"""

from . import scalespace
#from visuals import *
from scipy import ndimage
from scipy.ndimage import morphology
from skimage.feature import peak_local_max
import skimage

from matplotlib.pyplot import *


def blob_enhancement(im, blob_size):
    """Returns an image, where blob-sized blobs are enhanced. 
    Blobs are assumed to be darker than their surroundings here.
    """
    sigma = blob_size/5
    
    si, si_c = scalespace.shape_index(im, sigma)
    
    disk        = skimage.morphology.disk(0.5*blob_size)
    small_disk  = skimage.morphology.disk(0.25*blob_size)
    Y = contrast_enhancement(-intensity_stretch(si,3), disk)
    Y = contrast_enhancement(Y, disk)  # Repeat
    Y = morphology.filters.convolve(Y,small_disk/np.sum(small_disk)) # Circular average
    
    return Y, si
    
def local_minima(im, blob_radius):
    """Find local minima separated according to the given blob radius"""
    
    X = ndimage.filters.gaussian_filter(im, blob_radius/3)   # Gaussian smoothing
    
    extrema = peak_local_max(X, blob_radius)

    return extrema
    
def crop(im, xmin, xmax, ymin, ymax):
    return im[ymin:ymax+1,xmin:xmax+1]
    

def contrast_enhancement(im, se):
    """Contrast enhancement by 
    adding the original image to the top-hat filtered image, 
    and then subtract the bottom-hat filtered image.
    A disk shaped element with radius disk_radius is used for structuring element.
    Top hat filtering == white top hat filtering."""
    
    return im + morphology.white_tophat(im,footprint=se) - morphology.black_tophat(im,footprint=se)


def intensity_stretch(im, nstd):
    """Stretch intensities of image linearly between mean +/- nstd standard 
    deviations such that 1 = mu + nstd*sigma """
    
    mu = np.mean(im)
    sigma = np.std(im)
    old_min = mu - nstd*sigma
    old_max = mu + nstd*sigma
    im_out = (im - old_min)/(old_max-old_min)
    im_out[im<0] = 0
    im_out[im>1] = 1
    return im_out
    
    
def fitquadratic(im,center,n):
    """
    Fit a quadratic around a point and find its local minimum/maximum
    
    Center should be in row,column format
    """
    
    w = int((n-1)/2)
    start = center - w
    patch = im[start[0]:(start[0]+n), start[1]:(start[1]+n)]
    
    x0 = (n+1)/2
    x,y = np.meshgrid(range(n),range(n))
    x = x.flatten()+1-x0
    y = y.flatten()+1-x0
    
    # Constants for coefficient vector [f a b c d e]
    A = np.vstack([np.ones(n*n), x**2, y**2, x*y, x, y]).T
    coef = np.linalg.lstsq(A, patch.flatten())[0]   

    # Find minimum coordinates
    B = np.array([[2*coef[1], coef[3]], [coef[3], 2*coef[2]]])
    rhs = np.array([[-coef[4]],[-coef[5]]])  # [2a c; c 2b][x;y]=-[d;e]
    xymin = np.linalg.lstsq(B,rhs)[0]
    xymin += np.array([[center[1]],[center[0]]])
    return xymin, coef