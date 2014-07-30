# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 18:57:51 2014

@author: schackv
"""

from . import dm3lib_v099 as dm3lib

class dm3image():
    """This class wraps the dm3lib.DM3 class"""
        
    def __init__(self, image_path):
        """Initialize the DM3 image from the given path."""
        self.dm3f = dm3lib.DM3(image_path)
        
    def image(self):
        return self.dm3f.getImageData()
        
    def pixelsize(self):
        pixelsize, unit = self.dm3f.getPixelSize()
        return pixelsize, unit