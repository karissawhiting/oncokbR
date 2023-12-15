library(cbioportalR)
library(tidyverse)


# Get Data ----------------------------------------------------------

set_cbioportal_db("public")
all_gen <- get_genetics_by_study("blca_msk_tcga_2020")

mut_raw <- all_gen$mutation

mut <-  mut_raw %>%
  select("Hugo_Symbol" = hugoGeneSymbol,
        "Variant_Classification" = mutationType,
        "Chromosome" = chr,
        "Start_Position" = startPosition,
        "End_Position" = endPosition,
        "Reference_Allele" = referenceAllele) %>%
  mutate("Tumor_Seq_Allele2" =  NA)

cna <- all_gen$cna %>%
  rename(
    "Tumor_Sample_Barcode" = sampleId,
    "Hugo_Symbol" = hugoGeneSymbol,
         "Copy_Number_Alteration" = alteration)

sv <- all_gen$structural_variant %>%
  rename("GeneA" = site1HugoSymbol,
         "GeneB" = site2HugoSymbol)


write.table(mut, file = "/Users/kwhiting/Repositories/oncokb-annotator/data/blca_data/blca_mutation.txt", sep = '\t')
write.table(cna, file ="/Users/kwhiting/Repositories/oncokb-annotator/data/blca_data/blca_cna.txt", sep = '\t')
write.table(sv, file = "/Users/kwhiting/Repositories/oncokb-annotator/data/blca_data/blca_sv.txt", sep = '\t')

# Annotate Files ----------------------------------------------------------
library(tictoc)

#13.49 sec elapsed
tic()
mut_an <- annotate_mutations(mut_raw[1:100,])
toc()


#97.729 sec elapsed
tic()
mut_an <- annotate_mutations(mut_raw[1:1000,])
toc()

# 844.006 sec elapsed
tic()
mut_an <- annotate_mutations(mut_raw[1:10000,])
toc()

save(mut_an, file = here::here("data-raw", "blca_mutations_annotated.RData"))



mut_python_an <- read_table(file = "/Users/kwhiting/Repositories/oncokb-annotator/data/blca_data/blca_mutation.oncokb.txt")
write.table(cna, file = here::here("data-raw", "blca_cna.txt"), sep = '\t')
write.table(sv, file = here::here("data-raw", "blca_sv.txt"), sep = '\t')

# Read in Annotator Files -------------------------------------------------


