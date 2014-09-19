# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 15:43:43 2014

@author: jsve
"""

import os
import csv
import graphene.gridmatching

data_dir = r'E:\dtu\phd\graphene\simulated_images'


def run():
    files = _getfiles(os.path.join(data_dir,'simulations.csv'))
    
    for file in files:
        filepath = os.path.join(data_dir,file['Filename'])
        print(filepath)
        graphene.gridmatching.main(filepath)
        




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