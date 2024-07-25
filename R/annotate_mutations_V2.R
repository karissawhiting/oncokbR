library(purrr)
library(tidyverse)


hgvsg = c("7:g.140453136A>T", "3:g.179218294G>A")
referenceGenome = c("GRCh37", "GRCh37")

res <- annotate_mutations_by_hgvsg(hgvsg, referenceGenome)

resp_df<- map(res, ~bind_cols(.x))

