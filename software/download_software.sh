#!/bin/bash

cd software/
# Seurat V2
curl -L https://github.com/satijalab/seurat/archive/v2.3.3.tar.gz -o seurat.tar.gz
tar -xvzf seurat.tar.gz

# Seurat V3
git clone --single-branch --branch release/3.0 https://github.com/satijalab/seurat.git

git clone https://github.com/brianhie/scanorama.git
cd scanorama
git checkout b83e1ba87e6635825e155046b884dee2154c5d34
cd ..