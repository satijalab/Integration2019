import numpy as np
import sys
import os
from glob import glob
from argparse import ArgumentParser
from scipy.spatial import ConvexHull


parser = ArgumentParser(description='Find spatial positions of cells')
parser.add_argument('--starmap_path', help='path to starmap repository', required = True, type = str)
parser.add_argument('--data_path', help='path to root data folder', required = True, type = str)
options = parser.parse_args()


# From https://github.com/weallen/STARmap/blob/master/python/viz.py

def GetQHulls(labels):
    labels += 1
    Nlabels = labels.max()
    hulls = []
    coords = []
    num_cells = 0
    print('blah')
    for i in range(Nlabels):#enumerate(regionprops(labels)):
        print(i,"/",Nlabels)
        curr_coords = np.argwhere(labels==i)
        # size threshold of > 100 pixels and < 100000
        if curr_coords.shape[0] < 100000 and curr_coords.shape[0] > 1000:
            num_cells += 1
            hulls.append(ConvexHull(curr_coords))
            coords.append(curr_coords)
    print("Used %d / %d" % (num_cells, Nlabels))
    return hulls, coords

def get_polys(label_path, qhulls_file, centroids_file):
    labels = np.load(label_path)['labels']
    qhulls, coords = GetQHulls(labels)
    all_centroids = np.vstack([c.mean(0) for c in coords])
    np.savetxt(fname=centroids_file, X=all_centroids, delimiter="\t")
    with open(qhulls_file, 'w+') as outf:
        i = 1
        for l in qhulls:
            for x in l.vertices:
                outf.write(str(i) + "\t" + "\t".join(map(str, l.points[x].tolist())) + "\n")
            i += 1


def iterate_files(base_path, experiment):
    for mydir in os.listdir(base_path + experiment):
        abspath = (base_path + experiment + "/" + mydir)
        if os.path.isdir(abspath):
            get_polys(label_path = abspath + "/labels.npz",
                      qhulls_file = abspath + "/qhulls.tsv",
                      centroids_file = abspath + "/centroids.tsv")

iterate_files(base_path=options.data_path, experiment='visual_1020')
iterate_files(base_path=options.data_path, experiment='visual_160')
