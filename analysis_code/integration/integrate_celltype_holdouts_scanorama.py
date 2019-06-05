import sys, os
sys.path.insert(0, os.getcwd() + '/software/scanorama/bin')

from process import load_names, merge_datasets, save_datasets
from scanorama import *

from time import time
import numpy as np
import random

random.seed(42)
np.random.seed(42) 


NAMESPACE = 'panorama'

# Find the panoramas in the data.
def panorama(datasets_full, genes_list):
    if VERBOSE:
        print('Processing and reducing dimensionality...')
    datasets, genes = merge_datasets(datasets_full, genes_list)
    datasets_dimred, genes = process_data(datasets, genes)

    if VERBOSE:
        print('Finding panoramas...')
    panoramas = connect(datasets_dimred)

    if VERBOSE:
        print(panoramas)

    return panoramas

if __name__ == '__main__':
    from config import data_names
    
    # Load raw data from files.
    datasets, genes_list, n_cells = load_names(data_names)
    
    ##########################
    ## Panorama correction. ##
    ##########################

    # Put each of the datasets into a panorama.
    t0 = time()
    panoramas = panorama(datasets, genes_list)
    if VERBOSE:
        print('Found panoramas in {:.3f}s'.format(time() - t0))

    # Assemble and correct each panorama individually.
    t0 = time()
    for p_idx, pano in enumerate(panoramas):
        if VERBOSE:
            print('Building panorama {}'.format(p_idx))

        # Consider just the datasets in the panorama.
        pan_datasets = [ datasets[p] for p in pano ]
        pan_genes_list = [ genes_list[p] for p in pano ]

        # Do batch correction on those datasets.
        pan_datasets, genes = correct(pan_datasets, pan_genes_list)
        if VERBOSE:
            print('Found {} genes in panorama'.format(len(genes)))
            
        # Store batch corrected result.
        for i, p in enumerate(pano):
            datasets[p] = pan_datasets[i]
            genes_list[p] = genes
            
    if VERBOSE:
        print('Batch corrected panoramas in {:.3f}s'.format(time() - t0))

    save_datasets(datasets, genes, data_names)

