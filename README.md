# mitoClone

An R package for clonal tracking in single-cell RNA-seq data using nuclear and mitochondrial mutations. Based on Velten, Story et al., unpublished.

## Requirements and installation

This package has only been tested under Linux and Mac. It requires a current version of R (tested under 3.6.2) and python (tested under 3.7.3) and, importantly, an installation of gurobi and the gurobipy python package, [gurobi installation instructions](https://www.gurobi.com/documentation/9.0/quickstart_mac/software_installation_guid.html) and [instructions for installing gurobipy](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-]). gurobi is freely available for academic users,see http://www.gurobi.com 

If the use of conda is intended for managing the gurobi python environment, this conda environment is passed tp R using the parameter `python_env` of the `muta_cluster` function. See function documentation for detail, example: `python_env = "source activate gurobi"`.

Once the gurobipy package is available, installation of the package can be performed in less than one minute using
```
devtools::install_github("veltenlab/mitoClone")
```

## Demo and instructions

Demos and instructions are contained in the package vignettes that become available upon installation.

1. callingCohort: Instructions on how to filter mitochondrial mutations using the strategy applied in the manuscript (typical runtime: > 5 minutes on a single CPU)
2. calling:  Instructions on how to filter mitochondrial mutations if only data from a single individual is available (typical runtime: > 5 minutes on a single CPU)
3. clustering: Instructions on how to cluster mutations into a clonal hierarchy and how to assign cells to clones (typical runtime: > 2 hours on a single CPU)
4. denovo: Instructions on how to identify new mutations associated with the clones. This final part can have long runtimes depending on the dataset size and number of mutations. The vignette contains instructions that the user can adapt to the appropriate compute infrastructure.