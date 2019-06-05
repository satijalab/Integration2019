import sys, os
sys.path.insert(0, os.getcwd() + '/software/scanorama/bin')

import gzip
import numpy as np
from sklearn.preprocessing import normalize
from process import process_tab, process_mtx

if __name__ == '__main__':
    from config import data_names

    for name in data_names:
        if os.path.isdir(name):
            process_mtx(name, min_trans=0)
        elif os.path.isfile(name):
            process_tab(name, min_trans=0)
        elif os.path.isfile(name + '.txt'):
            process_tab(name + '.txt', min_trans=0)
        elif os.path.isfile(name + '.txt.gz'):
            process_tab(name + '.txt.gz', min_trans=0)
        else:
            sys.stderr.write('Warning: Could not find {}\n'.format(name))
        print('Successfully processed {}'.format(name))


