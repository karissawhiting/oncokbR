
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
[Chakravarty et
al. 2017](https://ascopubs.org/doi/full/10.1200/PO.17.00011).

{oncokbR} is a part of the
[gnomeverse](https://mskcc-epi-bio.github.io/genomeverse/), a collection
of R packages designed to work together seamlessly to create
reproducible clinico-genomic analysis pipelines.

## Installation

You can install the development version of oncokbR with:

``` r
remotes::install_github('karissawhiting/oncokbR')
```

## Authentication

In order to use this package, you will need to acquire an oncoKB API
token and save it to your .Renviron file (or wherever you store
credentials). You will need to register for an account on [OncoKB
Precision Oncology Knowledge Base](https://www.oncokb.org/). Once your
registration is approved, you will see a token under your [Account
Settings](https://www.oncokb.org/account/settings). Save it in your
.Renviron file as `ONCOKB_TOKEN`. This will save the token as an
environmental variable called `ONCOKB_TOKEN` that can be recognized in
all package functions.

*Tip: The following {usethis} function can easily find and open the
.Renviron for you:*

``` r
usethis::edit_r_environ()
#> ☐ Edit '/Users/kwhiting/.Renviron'.
#> ☐ Restart R for changes to take effect.
```

Paste the token you were given (using the format below) in the .Renviron
file and save the file changes. After saving you should restart your R
session to ensure the token is saved and recognized.

``` r
ONCOKB_TOKEN = 'YOUR_TOKEN'
```

## Annotate Data

### Annotate Mutations:

``` r
library(oncokbR)
library(dplyr)
#> Warning: package 'dplyr' was built under R version 4.2.3
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

Annotate MAF data with tumor type indicated for annotations on
oncogenicity and oncoKB treatment levels. This will show the highest
sensitivity levels for bladder cancer specifically:

``` r
blca_mutation <- oncokbR::blca_mutation %>%
  mutate(tumor_type = "BLCA")

annotated_tt <- annotate_mutations(mutations = blca_mutation[1:50,])
#> ✔ Token valid for "whitingk@mskcc.org"
#> ℹ Token expiration date: 2025-04-20

annotated_tt$oncokb_oncogenic 
#>  [1] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#>  [5] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#>  [9] "Unknown"          "Oncogenic"        "Unknown"          "Unknown"         
#> [13] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#> [17] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#> [21] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#> [25] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#> [29] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#> [33] "Unknown"          "Likely Neutral"   "Likely Neutral"   "Unknown"         
#> [37] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#> [41] "Likely Oncogenic" "Unknown"          "Unknown"          "Unknown"         
#> [45] "Unknown"          "Likely Oncogenic" "Unknown"          "Unknown"         
#> [49] "Unknown"          "Unknown"
```

``` r
annotated_tt %>%
  select(oncokb_highest_sensitive_level) %>% 
  table()
#> oncokb_highest_sensitive_level
#> LEVEL_1 
#>       1
```

You can also annotated with no tumor type data for oncogenicity. Here
you’ll see the highest sensitive lever is higher than in the bladder
specific annotation because it’s referring to a non bladder cancer type
that has a higher level.

``` r
blca_mutation <- oncokbR::blca_mutation

annotated_no_tt <- annotate_mutations(mutations = blca_mutation[1:50,])
#> ✔ Token valid for "whitingk@mskcc.org"
#> ℹ Token expiration date: 2025-04-20
#> ℹ No "tumor_type" found in data. No treatment-level annotations will be returned.

annotated_no_tt$oncokb_oncogenic
#>  [1] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#>  [5] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#>  [9] "Unknown"          "Oncogenic"        "Unknown"          "Unknown"         
#> [13] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#> [17] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#> [21] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#> [25] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#> [29] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#> [33] "Unknown"          "Likely Neutral"   "Likely Neutral"   "Unknown"         
#> [37] "Unknown"          "Unknown"          "Unknown"          "Unknown"         
#> [41] "Likely Oncogenic" "Unknown"          "Unknown"          "Unknown"         
#> [45] "Unknown"          "Likely Oncogenic" "Unknown"          "Unknown"         
#> [49] "Unknown"          "Unknown"
```

``` r
annotated_no_tt %>%
  select(oncokb_highest_sensitive_level) %>% 
  table()
#> oncokb_highest_sensitive_level
#> LEVEL_1 
#>       1
```

### Annotate CNA:

``` r
blca_cna <- blca_cna %>%
  mutate(tumor_type = "BLCA")


annotated <- annotate_cna(blca_cna[1:10,])
#> ✔ Token valid for "whitingk@mskcc.org"
#> ℹ Token expiration date: 2025-04-20
table(annotated$oncokb_oncogenic)
#> 
#> Likely Oncogenic        Oncogenic          Unknown 
#>                3                6                1
```

### Annotate Structural Variants:

``` r

blca_sv <- blca_sv %>%
  mutate(tumor_type = "BLCA")


annotated <- annotate_sv(blca_sv[1:10,])
#> ✔ Token valid for "whitingk@mskcc.org"
#> ℹ Token expiration date: 2025-04-20
table(annotated$oncokb_oncogenic)
#> 
#> Likely Oncogenic          Unknown 
#>                5                5
```

## Contributing

Please note that {oncokbR} is released with a [Contributor Code of
Conduct](https://www.karissawhiting.com/oncokbR/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

Thank you to contributors!

[@david-nemirovsky](https://github.com/david-nemirovsky),
[@patrick120161](https://github.com/patrick120161),
[@edrill](https://github.com/edrill),
[@jalavery](https://github.com/jalavery),
[@stl2137](https://github.com/stl2137)
