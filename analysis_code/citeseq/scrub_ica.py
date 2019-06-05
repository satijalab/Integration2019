#! /usr/bin/env python3 

import scrublet as scr
import scipy.io
import numpy as np 


def run_scrublet(infile):
    """Run scrublet on a count matrix
    """
    print('Processsing '+infile)
    E = scipy.io.mmread(infile)
    E = E.T.tocsc()

    # filtering/preprocessing parameters:
    min_counts = 1
    min_cells = 1
    vscore_percentile = 85
    n_pc = 30

    # doublet detector parameters:
    expected_doublet_rate = 0.06 
    sim_doublet_ratio = 3
    n_neighbors = 50

    scrublet_results = scr.compute_doublet_scores(
        E, 
        min_counts = min_counts, 
        min_cells = min_cells, 
        vscore_percentile = vscore_percentile, 
        n_prin_comps = n_pc,
        scaling_method = 'zscore',
        expected_doublet_rate = expected_doublet_rate,
        sim_doublet_ratio = sim_doublet_ratio,
        n_neighbors = n_neighbors, 
        use_approx_neighbors = True, 
        get_doublet_neighbor_parents = False
    )

    outfile = infile.split('.mtx')[0]
    with open(outfile+'_doublets.tsv', "w+") as outf:
        obs = scrublet_results['doublet_scores_observed_cells'].tolist()
        sim = scrublet_results['doublet_scores_simulated_doublets'].tolist()
        for i in range(len(obs)):
            outf.write(str(obs[i]) + "\t" + str(sim[i]) + "\n")

if __name__ == '__main__':
    import glob
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Run scrublet given a path to files ending in .mtx')
    parser.add_argument('--path', help='Path to count matrix files', required = True)
    options = parser.parse_args()
    files = glob.glob(options.path+"/*.mtx")
    for i in files:
        run_scrublet(infile = i)
    
