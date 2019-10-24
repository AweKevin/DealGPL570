
# DealGPL570

<!-- badges: start -->
<!-- badges: end -->

The goal of DealGPL570 is to deal GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array) RAW.tar file using the robust multi-array average expression measure.
#'

## Installation

You can install the released version of DealGPL570 from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("DealGPL570")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(DealGPL570)
file <- system.file("extdata", "GSE104683_RAW.tar", package = "DealGPL570")
result <- DealGPL570(file = file)
```

