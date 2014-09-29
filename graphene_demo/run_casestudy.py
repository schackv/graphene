# -*- coding: utf-8 -*-
"""
Processes all exposures of the two cases 'hole' and 'regular' (or 'pristine')
and creates CDF plots of these fits.

To run this, you need to point to two folders with .dm3 files (or other image files).
If you want to run for just a single folder, just comment out the part where the 
comparing CDF plots are made.

Created on Mon Sep 29 07:32:16 2014

@author: jsve
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 15:43:43 2014

@author: jsve
"""

import os
import glob
import logging
import graphene.gridmatching
from graphene import misc
import matplotlib.pyplot as plt

data_dirs = ['./hole_data','./regular_data']    # Suggestion for use: Create symbolic links with these names
settings_files = ['hole_settings.txt','regular_settings.txt']

# Field names from csv file
NM_X_LABEL = 'Horizontal resolution [nm/px]'
NM_Y_LABEL = 'Vertical resolution [nm/px]'
VERTICAL_LABEL = 'Vertical'
OTHER_LABEL = 'Other'

def run(data_dir,settings_file,do_plot=True):
    files = _getfiles(data_dir)
    
    # Run analysis of all exposures
    bond_lengths = []
    for filepath in files:
        out = filepath + '_results'
        
        graphene.gridmatching.main(filepath, output_dir=out, force_fresh=False, settings_file = settings_file, do_plot=False)
        
        # Read result and get stats
        info = misc._readdict(os.path.join(out, 'info.txt'))
        H = graphene.grid.HexagonalGrid.from_textfile( os.path.join(out, 'atoms') )
#        bins, oriented_lengths = H.edge_lengths_by_orientation()
        bond_lengths.append(H.edge_lengths(scale_xy=info['nm_per_pixel_est']))
        
        logging.info('Processed {:s}'.format(filepath))
    
    if do_plot:
        # Make distribution plots of all exposures' bond lengths  
        plt.figure(figsize=(10,6))
        for L in bond_lengths:
            xcdf, F = misc.ecdf(L)
            plt.plot(xcdf,F)
        plt.xlim((0.117, 0.167))
        plt.grid()
        plt.savefig(os.path.join(data_dir,'cdf.png'))
        #    plt.show(block=True)
    
    return bond_lengths

def _getfiles(folder,pattern='*.dm3'):
    files = glob.glob(os.path.join(folder, pattern))
    return files

def cdf_plot(bondlength_dict,colors):
    """Make CDF plot of all exposures for the two cases.
    """
    plt.figure()
    for key, lengths in bondlength_dict.items():
        for i, L in enumerate(lengths):
            xcdf, F = misc.ecdf(L)
            if i==1:
                plt.plot(xcdf,F,color=colors[key],label=key)
            else: 
                plt.plot(xcdf,F,color=colors[key])
        
    plt.legend()
    plt.xlim((0.117, 0.167))
    plt.grid()
    plt.show(block=True)
        
        

if __name__=='__main__':
    hole = run(data_dirs[0],settings_files[0])     # Run for 'hole'
    regular = run(data_dirs[1],settings_files[1])     # Run for 'regular'
    
    cdf_plot({'Hole': hole, 'Pristine': regular},colors={'Hole': 'r', 'Pristine': 'b'})
    
    
    
    
    
    
    
    
    
    
    
    