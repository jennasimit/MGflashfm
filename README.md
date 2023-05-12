

<!-- README.md is generated from README.Rmd. Please edit that file -->



# MGflashfm

<!-- badges: start -->
<!-- badges: end -->

The goal of MGflashfm is to use GWAS summary statistics to jointly fine-map genetic associations for several 
related quantitative traits across multiple population groups. MGfm is also available to fine-map genetic associations
for a single quantitative trait across multiple population groups.

Website available at: https://jennasimit.github.io/MGflashfm/

We have applied these methods to GWAS results from four lipids traits and five population groups, as made available by
 the Global Lipids Genetic Consortium (GLGC). Our analysis scripts are available here:

https://github.com/fz-cambridge/MGflashfm-GLGC-analysis

## System Requirements

MGflashfm could be installed with ease on versions of R > 4.2.1.
If installing on a Windows machine, Rtools must be installed.
Installation time is estimated as 2 minutes.

## Installation Guide

## Short version

``` r
# install.packages("devtools")
devtools::install_github("jennasimit/MGflashfm")
```

## Longer version (if above fails)

The following packages from CRAN and Bioconductor are required:

``` r
install.packages("parallel")
install.packages("Matrix")
install.packages("gtools")
install.packages("rlist")
```

as well as the following dependencies from GitHUb


``` r
# install and load flashfm and R2BGLiMS
remotes::install_github("jennasimit/flashfm")
remotes::install_github("pjnewcombe/R2BGLiMS")
```

NB: Must have a Java JDK installed in order to install and run R2BGLiMS. This is only needed if you need to run single-trait fine-mapping using JAM. 
If single-trait fine-mapping results are available, then it is not necessary to have Java JDK installed.

``` r
remotes::install_github("jennasimit/MGflashfm")
library(MGflashfm)
library(R2BGLiMS)  # if running internal JAM functions for single-trait fine-mapping
```


