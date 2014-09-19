# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 14:53:49 2014

@author: jsve
"""

import os
from graphene import gridmatching
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f","--file",help="Image filename",required=True)
parser.add_argument("-s","--settings_file",help="Parameter settings file. If not specified, default values will be used.")
parser.add_argument("-nm","--nm_per_pixel", default=(0.01153, 0.01553), type=float,nargs=2, metavar=('x','y'), help="Image resolution in nanometer per pixel (default: %(default)s). Use -nm -1 -1 to estimate resolution from image by assuming known CC bond length.")
parser.add_argument("-o","--output_dir",help="Output directory (defaults to FILE_results/")
parser.add_argument("-ff","--force-fresh",dest='force_fresh',action='store_true',help='Force processing of file, even though results already exist')
parser.add_argument("-p","--plot", dest='do_plot',action='store_true',help="Enable output of plots.")
parser.add_argument("-np","--no-plot", dest='do_plot',action='store_false',help="Disable output of plots.")
parser.add_argument("-pd","--plot_dir",help="Plot output directory (defaults to OUTPUT_DIR/plots/). Only relevant when --plot is used.")
parser.add_argument("-l","--loglevel",default="INFO",help="Logging level (default: %(default)s)",choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'])
parser.set_defaults(do_plot=False,force_fresh=False)
parser.print_help()
args=parser.parse_args()

if (args.output_dir==None):
    args.output_dir = '{:s}_results'.format(args.file)
if (args.plot_dir==None):
    args.plot_dir = os.path.join(args.output_dir,'plots')
if any( nm < 0 for nm in args.nm_per_pixel):
    args.nm_per_pixel = None

if __name__ == '__main__':
    gridmatching.main(**vars(args))