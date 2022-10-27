
<!-- README.md is generated from README.Rmd. Please edit that file -->

# oncokbR

<!-- badges: start -->
<!-- badges: end -->

{oncokbR} allows you to annotate genomic data sets using MSK’s [OncoKB
Precision Oncology Knowledge Base](https://www.oncokb.org/) directly in
R. The package wraps existing [OncoKB API
endpoints](https://api.oncokb.org/oncokb-website/api) so R users can
easily leverage the existing API to annotate data on mutations, copy
number alterations and fusions.

This package is compatible with oncoKB data v3.16, but is subject to
change as new versions are released. For more information on oncoKB, see
\[Chakravarty et
al. 2017\[(<https://ascopubs.org/doi/full/10.1200/PO.17.00011>).

Gao et al. Sci. Signal. 2013 Cerami et al. Cancer Discov. 2012

## Installation

You can install the development version of oncokbR like so:

``` r
remotes::install_github('karissawhiting/oncokbR")
```

## Annotate Data

Annotate Mutations:

``` r
library(oncokbR)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

blca_mutation <- oncokbR::blca_mutation %>%
  mutate(tumor_type = "ACC")


annotated <- annotate_mutations(blca_mutation[1:20,])
table(annotated$oncogenic)
#> 
#> Oncogenic   Unknown 
#>         1        19
```

Annotate CNA:

``` r
blca_cna <- blca_cna %>%
  mutate(tumor_type = "ACC")


annotated <- annotate_cna(blca_cna[1:10,])
table(annotated$oncogenic)
#> 
#> Likely Oncogenic        Oncogenic          Unknown 
#>                3                6                1
```