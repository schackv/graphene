# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 15:43:09 2014

@author: jsve
"""


import os
import shutil
from copy import deepcopy
import logging
import datetime
import matplotlib.pyplot as plt
import numpy as np
from . import options, misc, grid, bgm, imtools, lattice, alternating_graphcut


INFO_FILENAME = r'info.txt'

def _saveplot(name,formats={'png','pdf'}):
    for f in formats:
        fullname = name + '.' + f
        logging.debug('Writing {:s}'.format(fullname))
        plt.savefig(fullname)

def main(file, settings_file=None, nm_per_pixel=None,output_dir=None, force_fresh=False, do_plot=False, plot_dir=None, loglevel="INFO"):
    
    # Set logging level
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: {:s}'.format(loglevel))
    logging.basicConfig(level=numeric_level)
    
    os.makedirs(output_dir,exist_ok=True)
    # Read settings file and write immediately to output dir
    logging.info('Reading settings')
    config = _readconfig(settings_file)
    misc._writedict(os.path.join(output_dir,'options.txt'),config)
    
    # Process (and save output)
    if not os.path.exists(os.path.join(output_dir,INFO_FILENAME)) or force_fresh:
        logging.info('Processing {:s}'.format(file))
        process_image(file,output_dir=output_dir,opts=config)
        logging.info('Done.')
    else:
        logging.info('Skipping processing of {:s} - results already exist.'.format(file))
    
    # Make plots (if chosen)
    if do_plot:
        logging.info('Saving plots to {:s}'.format(plot_dir))
        make_plots(output_dir,plot_dir, nm_per_pixel=nm_per_pixel)
        logging.info('Done.')
    
    
def process_image(filename, output_dir=None, opts={}):
    
    # Read image
    im = imtools.read_image(filename)
    lp = lattice.parameters()

    lp.compute(im,opts['lattice'])
    logging.debug('Estimated hexagonal side length in pixels = {:8f}'.format(lp.t))
    
    ## Initial points
    extrema, _, _= lattice.hexagonal_centers(im, lp.t)
    logging.debug('{:g} initial points detected.'.format(len(extrema)))
    
    ## Initial grid
    xy, simplices = alternating_graphcut.cleanup(extrema,im)
    logging.debug('Alternating graph cut keeps {:g} points.'.format(len(xy)))
    G_initial = grid.TriangularGrid.from_simplices(xy,simplices)
    
    ## Fine adjustment
    logging.debug('Now fine adjusting...')
    G = deepcopy(G_initial)
    cs = bgm.welsh_powell(G.edges())
    model = bgm.AdaptiveGrid(im, G.edges(), G.simplices) 
    xy_hat, E, history = bgm.fit_grid(im,G_initial.xy, G_initial.edges(), coding_scheme=cs, beta=opts['fine_adjustment']['beta'], gridenergydefinition=model, anneal_opts=opts['annealing'])
    G.xy = xy_hat
    logging.debug('Fine adjustment complete.')

    ## Place atoms
    H = grid.HexagonalGrid.from_triangular(G)
    logging.debug('{:g} atoms placed.'.format(len(H.xy)))
    
    ## Save data and output path
    if output_dir != None:
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        filename_base = os.path.basename(filename)
        shutil.copyfile(filename, os.path.join(output_dir,filename_base)) # Copy image
        
        # Write grids to txt files
        G_initial.write(os.path.join(output_dir,'tri_initial'))
        G.write(os.path.join(output_dir,'tri_final'))
        H.write(os.path.join(output_dir,'atoms'))
        
        # Save info-file
        info = {'filename': filename_base,'timestamp': str(datetime.datetime.now()), 'nm_per_pixel_est': lp.nm_per_pixel_est()}
        misc._writedict(os.path.join(output_dir,INFO_FILENAME),info)
        
        logging.info('Results written to {:s}'.format(output_dir))
        

def make_plots(result_dir,plot_dir,nm_per_pixel=None):
    os.makedirs(plot_dir,exist_ok=True)
    
    info = misc._readdict(os.path.join(result_dir,INFO_FILENAME))
    
    if nm_per_pixel is None:
        nm_per_pixel = info['nm_per_pixel_est']
        logging.debug('Resolution not supplied. Using estimate of {:.6f} nm per pixel.'.format(nm_per_pixel))
    
    pname = lambda s: os.path.join(plot_dir,s)
    
    # Read image
    im = imtools.read_image(os.path.join(result_dir,info['filename']))
    
    # Write image of fitted grid
    G = grid.Grid.from_textfile(os.path.join(result_dir,'tri_final'))
    plt.figure()
    plt.imshow(im,cmap=plt.cm.gray,interpolation='nearest')
    G.plot()
    _saveplot(pname('tri_final'))
    plt.close()
    
    # Write hexagonal grid
    H = grid.HexagonalGrid.from_textfile(os.path.join(result_dir,'atoms'))
    plt.figure()
    plt.imshow(im,cmap=plt.cm.gray,interpolation='nearest')
    H.plot()
    _saveplot(pname('atoms'))
    plt.close()
    
    # CDF of bond lengths
    lengths_nm = H.edge_lengths(nm_per_pixel)
    plt.figure()
    cdf_plot(lengths_nm)
    _saveplot(pname('cdf_nm'))
    plt.close()
    
    # CDF of oriented bond lengths
    orientations = H.edge_orientations()
    bin_centers, bin_idx = misc.circular_binning(orientations,half_circle=True)
    for b, bin_center in enumerate(bin_centers):
        oriented_lengths = [x for (x,idx) in zip(lengths_nm,bin_idx) if idx==b]
        plt.figure()
        cdf_plot(oriented_lengths)
        bin_str = '{:.0f}'.format(np.rad2deg(bin_center))
        plt.title(bin_str)
        _saveplot(pname('o{:s}_cdf_nm'.format(bin_str)))
        plt.close()
        print('Orientation {:.2f}: {:.5f} +/- {:.5f}'.format(np.rad2deg(bin_center),np.mean(oriented_lengths), np.std(oriented_lengths)))

def cdf_plot(vals):
    xcdf, F = misc.ecdf(vals)
    plt.plot(xcdf,F)
    plt.xlim((0.127, 0.157))
    plt.grid()


## Config file handling
def _readconfig(filename):
    """Returns a dictionary of the configuration."""
    # Default options
    config = options.defaults
    # Merge with other options file
    if filename != None:
        new_options = misc._readdict(filename)
        # Union default and new options
        config = misc._merge(config,new_options)
    
    logging.debug(config)
    
    return config




#    # Profiling
#    import cProfile
#    cProfile.run('main(**vars(args))','stats')
#    
#    import pstats
#    p = pstats.Stats('stats')
#    p.strip_dirs().sort_stats(-1).print_stats()
#    
#    p.sort_stats('cumulative').print_stats(20)
    
    
    