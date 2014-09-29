# -*- coding: utf-8 -*-
"""
Simulation studies processing. Fits a grid to each image, extracts bond lengths,
calculates deviation from expected bond length and prints out a table with results.

The expected bond lengths expect that the simulated images are scaled along the 
y-axis only, but this can easily be modified.

NB: Depends on simulated images and a .csv file with the following column names
for each simulated image:

'Filename':                         Filename of image
'GroupId':                          Integer, only used for sorting the final table (low to high)
'Horizontal resolution [nm/px]':    Image resolution in x-direction
'Vertical resolution [nm/px]':      Image resolution in y-direction
'Vertical':                         Expected bond length for the vertical bonds
'Other':                            Expected bond lengths for the two other bond orientations

Created on Fri Sep 19 15:43:43 2014
@author: jsve
"""

import os
import csv
import numpy as np
import graphene.gridmatching
from tabulate import tabulate

data_dir = './simulated_images'     # Suggestion for use: Create a symbolic link with this name to the actual folder

# Field names from csv file
NM_X_LABEL = 'Horizontal resolution [nm/px]'
NM_Y_LABEL = 'Vertical resolution [nm/px]'
VERTICAL_LABEL = 'Vertical'
OTHER_LABEL = 'Other'

def run():
    files = _getfiles(os.path.join(data_dir,'simulations.csv'))

    # Run analysis of all simulated images
    table = []
    for file in sorted(files,key=lambda f: int(f['GroupId'])):
        filepath = os.path.join(data_dir,file['Filename'])
        out = filepath + '_results'
        
        graphene.gridmatching.main(filepath, output_dir=out,force_fresh=False)
        
        # Read result and get stats
        nm_per_pixel = float(file[NM_X_LABEL]), float(file[NM_Y_LABEL])
        H = graphene.grid.HexagonalGrid.from_textfile( os.path.join(out, 'atoms') )
        bins, oriented_lengths = H.edge_lengths_by_orientation(scale_xy=nm_per_pixel)
        
        # Get expected lengths
        expected_lengths = [float(file[VERTICAL_LABEL]) if b==90 else float(file[OTHER_LABEL]) for b in np.rad2deg(bins)]
        
        row = [file['Filename']]
        row.extend(['{:.4f} +/- {:5f}'.format(np.mean(l-e),np.std(l-e)) for l, e in zip(oriented_lengths,expected_lengths)])
        table.append(row)
        
    print(tabulate(table))
        

def _getfiles(specfile):
    files = []
    with open(specfile,'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=';')
        columns = csvreader.__next__()
        for row in csvreader:
            file = dict(zip(columns,row))
            files.append(file)
    return files

if __name__=='__main__':
    run()