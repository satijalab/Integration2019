# Comprehensive integration of single-cell data

This repository provides code used to perform analyses presented in [*Comprehensive integration of single-cell data*](https://doi.org/10.1016/j.cell.2019.05.031).

An open-access bioRxiv version of the manuscript can be found [here](https://doi.org/10.1101/460147).

## How to run the analysis

All major components of the analysis can be reproduced using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workfow provided in this repository. Please note that there is a substantial amount of data required to reproduce the entire analysis (over 100 Gb).

Individual plots can be reproduced by setting the snakemake target. For example, figure 2 of the paper can be reproduced by running:

```
snakemake figure2
```

All the available rules can be viewed in the snakemake file.

## Citation

```
@article{,
  title    = "Comprehensive integration of single-cell data",
  author   = "Stuart, Tim and Butler, Andrew and Hoffman, Paul and Hafemeister,
              Christoph and Papalexi, Efthymia and Mauck, William M and Hao, Yuhan and
              Stoeckius, Marlon and Smibert, Peter and Satija, Rahul",
  journal  = "Cell",
  month    =  jun,
  year     =  2019,
  url      =  "https://doi.org/10.1016/j.cell.2019.05.031"
  language = "en"
}
```