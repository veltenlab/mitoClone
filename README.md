# mitoClone R Package

The tool is used for performing the analysis of clonal heterogeneity based on nuclear and mitochondrial mutations in single cell RNA or DNA sequencing.

## 1. System Requirements:
   - Linux/Mac
   - R 3.5+
   - Python 2.7, 3.6, or 3.7
   - Gurobi 9.0.0+
   
Importantly, an installation of both Gurobi and the gurobipy python package, [Gurobi Installation Instructions](https://www.gurobi.com/documentation/9.0/quickstart_mac/software_installation_guid.html) and [Instructions for Installing gurobipy](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-]). Gurobi is freely available for academic users see http://www.gurobi.com.

See **DESCRIPTION** file for specific R package requirements.

The software has been successfully implemented and tested using: Python 3.6.5, R 4.0.0, and Gurobi 9.0.3 on CentOS 7.

## 2. Installation
For manual package installation use the command:

`git clone https://github.com/veltenlab/mitoClone.git`

For installing the library from github into R directly use:

``` r
library(devtools)
devtools::install_github("veltenlab/mitoClone", build_vignettes = FALSE)
```

Estimated installation time: < 1 hour*

*\*Not including acquiring and installing Gurobi license.*

## 3. Demo

Please see R vignettes for further instructions and a demo using real data. Use the command `vignette("mitoClone")` after loading the library (see Instructions) to list all available tutorials.

Following included tutorials should output figures reproducing publication results.

Estimated demo completion time: < 1 hour

## 4. Usage Instructions

After installing all dependencies, open an R session and load the library from its installation directory (e.g. mitoClone-master) using the following command:

``` r
library(devtools)
load_all('mitoClone-master')
```

Or if installed via `install_github`:

``` r
library(mitoClone)
```

Please make sure to set your environmental python variables correctly for use of gurobi. See for example the `python_env` parameter used in the `muta_cluster` function.

Again please view the R vignettes for usage possibilities. See the following webpages (located in the cloned github folder `inst/doc`) for various tutorials.

   - **callingCohort**: Instructions on how to filter mitochondrial mutations using the strategy applied in the manuscript (typical runtime: > 5 minutes)
   - **calling**: Instructions on how to filter mitochondrial mutations if only data from a single individual is available (typical runtime: > 5 minutes )
   - **clustering**: Instructions on how to cluster mutations into a clonal hierarchy and how to assign cells to clones (typical runtime: > 2 hours)
   - **denovo**: Instructions on how to identify new mutations associated with the clones. This final part can have long runtimes depending on the dataset size and number of mutations. The vignettes contain further instructions so that user can adapt parameters to their appropriate computational infrastructure.

