# mitoClone R Package

The tool is used for performing the analysis of clonal heterogeneity based on nuclear and mitochondrial mutations in single cell RNA or DNA sequencing.

## 1. System requirements:
   - Linux/Mac OS
   - R 3.5+
   - Python 2.7, 3.6, or 3.7
   - Gurobi 9.0.0+

See DESCRIPTION file for specific R package requirements.

The software has been successfully implemented  and tested successfully using: Python 3.6.5, R 4.0.0 , and Gurobi 9.0.3 on CentOS 7.

## 2. Installation
For manual package installation use the command:

`git clone https://github.com/veltenlab/mitoClone.git`

For installing the library from github into R directly use:

``` r
library(devtools)
devtools::install_github("veltenlab/mitoClone", build_vignettes = TRUE)
```

Estimated installation time: < 1 hour*

*Not including acquiring and installing Gurobi license.

## 3. Demo

Please see R vignettes for demo and further instructions. Use the command vignette("mitoseq") after loading the library (see Instructions) to list all available tutorials.

Estimated demo completion time: < 1 hour

## 4. Usage Instructions

After installing all dependencies, open an R session and load the library from its installation directory (e.g. mitoClone-master) using the following command:

``` r
library(devtools)
load_all(mitoClone-master)
```

Or if installed via `install_github`:

``` r
library(mitoseq)
```

Please make sure to set your environmental python variables correctly for use of gurobi. See the `python_env` parameter used in the `muta_cluster` function.

Again please view the R vignettes for usage possibilities. See the following webpages (located in the cloned github folder) for various tutorials.


**callingCohort**: Instructions on how to filter mitochondrial mutations using the strategy applied in the manuscript (typical runtime: > 5 minutes on a single CPU)

**calling**: Instructions on how to filter mitochondrial mutations if only data from a single individual is available (typical runtime: > 5 minutes on a single CPU)

**clustering**: Instructions on how to cluster mutations into a clonal hierarchy and how to assign cells to clones (typical runtime: > 2 hours on a single CPU)

**denovo**: Instructions on how to identify new mutations associated with the clones. This final part can have long runtimes depending on the dataset size and number of mutations. The vignette contains instructions that the user can adapt to the appropriate compute infrastructure.

