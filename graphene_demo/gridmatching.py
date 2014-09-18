# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 15:43:09 2014

@author: jsve
"""


import os
import shutil
from copy import deepcopy
import argparse
import logging
import json
import datetime
from graphene import *
import matplotlib.pyplot as plt




INFO_FILENAME = r'info.txt'

def _saveplot(name,formats={'png'}):
    for f in formats:
        fullname = name + '.' + f
        logging.debug('Writing {:s}'.format(fullname))
        plt.savefig(fullname)

def main(file, settings_file=None, output_dir=None, force_fresh=False, do_plot=False, plot_dir=None, loglevel="INFO"):
    
    # Set logging level
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: {:s}'.format(loglevel))
    logging.basicConfig(level=numeric_level)
    
    os.makedirs(output_dir,exist_ok=True)
    # Read settings file and write immediately to output dir
    logging.info('Reading settings')
    config = _readconfig(settings_file)
    _writedict(os.path.join(output_dir,'options.txt'),config)
    
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
        make_plots(output_dir,plot_dir)
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
        info = {'filename': filename_base,'timestamp': str(datetime.datetime.now())}
        _writedict(os.path.join(output_dir,INFO_FILENAME),info)
        
        logging.info('Results written to {:s}'.format(output_dir))
        

def make_plots(result_dir,plot_dir,nm_per_pixel=0.0115):
    os.makedirs(plot_dir,exist_ok=True)
    
    info = _readdict(os.path.join(result_dir,INFO_FILENAME))
    
    pname = lambda s: os.path.join(plot_dir,s)
    
    # Read image
    im = imtools.read_image(os.path.join(result_dir,info['filename']))
    
    # Write image of fitted grid
    G = grid.Grid.from_textfile(os.path.join(result_dir,'tri_final'))
    plt.figure()
    plt.imshow(im,cmap=plt.cm.gray,interpolation='nearest')
    G.plot()
    _saveplot(pname('tri_final'))
    
    # Write hexagonal grid
    H = grid.HexagonalGrid.from_textfile(os.path.join(result_dir,'atoms'))
    plt.figure()
    plt.imshow(im,cmap=plt.cm.gray,interpolation='nearest')
    H.plot()
    _saveplot(pname('atoms'))
    
    # CDF of bond lengths
    lengths_px = H.edge_lengths()
    lengths_nm = lengths_px * nm_per_pixel
    plt.figure()
    xcdf, F = misc.ecdf(lengths_nm)
    plt.plot(xcdf,F)
    plt.show(block=True)
    _saveplot(pname('cdf_nm'))
            
## Config file handling
def _readconfig(filename):
    """Returns a dictionary of the configuration."""
    # Default options
    config = options.defaults
    # Merge with other options file
    if filename != None:
        new_options = _readdict(filename)
        # Union default and new options
        config = _merge(config,new_options)
    
    logging.debug(config)
    
    return config
    
def _readdict(filename):
    with open(filename, 'r') as f:
        d = json.load(f)
    return d

def _writedict(filename,config_dict):
    """Write a dictionary to file."""
    
    with open(filename,'w') as f:
        json.dump(config_dict,f,indent=4)
        
def _merge(a,b,path=None):
    "merges b into a"
    if path is None: path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                _merge(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass # same leaf value
            else:   # a values are overwritten by b values
                a[key] = b[key]
#                raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file",help="Image filename",required=True)
    parser.add_argument("-s","--settings_file",help="Parameter settings file. If not specified, default values will be used.")
    parser.add_argument("-o","--output_dir",help="Output directory (defaults to FILE_results/")
    parser.add_argument("-ff","--force-fresh",dest='force_fresh',action='store_true',help='Force processing of file, even though results already exist')
    parser.add_argument("-p","--plot", dest='do_plot',action='store_true',help="Enable output of plots.")
    parser.add_argument("-np","--no-plot", dest='do_plot',action='store_false',help="Disable output of plots.")
    parser.add_argument("-pd","--plot_dir",help="Plot output directory (defaults to OUTPUT_DIR/plots/). Only relevant when --plot is used.")
    parser.add_argument("-l","--loglevel",default="INFO",help="Logging level (default: %(default)s)",choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'])
    parser.set_defaults(do_plot=False,force_fresh=False)
#    parser.print_help()
    args=parser.parse_args()
    
    if (args.output_dir==None):
        args.output_dir = '{:s}_results'.format(args.file)
    if (args.plot_dir==None):
        args.plot_dir = os.path.join(args.output_dir,'plots')
    
    main(**vars(args))
#    
#    # Profiling
#    import cProfile
#    cProfile.run('main(**vars(args))','stats')
#    
#    import pstats
#    p = pstats.Stats('stats')
#    p.strip_dirs().sort_stats(-1).print_stats()
#    
#    p.sort_stats('cumulative').print_stats(20)
    
    
    
