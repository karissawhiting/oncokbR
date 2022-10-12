library(cbioportalR)

set_cbioportal_db("public")
blca_genetics <- get_genetics_by_study("blca_nmibc_2017")

blca_mutation <- blca_genetics$mutation
blca_cna <- blca_genetics$cna
blca_sv <- blca_genetics$structural_variant

usethis::use_data(blca_mutation, overwrite = TRUE)
usethis::use_data(blca_cna, overwrite = TRUE)
usethis::use_data(blca_sv, overwrite = TRUE)


set_cbioportal_db("public")
acc_genetics <- get_genetics_by_study("acc_tcga")

acc_mutation <- acc_genetics$mutation
acc_cna <- acc_genetics$cna
#acc_sv <- acc_genetics$structural_variant

usethis::use_data(acc_mutation, overwrite = TRUE)
usethis::use_data(acc_cna, overwrite = TRUE)

