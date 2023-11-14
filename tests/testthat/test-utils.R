

test_that("select column returns work", {
  simp <- annotate_mutations(oncokbR::blca_mutation[1:10, ],
                             return_simple = TRUE,
                             return_query_params = FALSE)

  should_be_in_res <- output_dictionary$output_column_name[
    which(output_dictionary$include_in_simple_output == "yes")
  ]

  expect_equal(names(select(simp, contains("oncokb_"))), should_be_in_res)

  simp_query <- annotate_mutations(oncokbR::blca_mutation[1:10, ],
                                   return_simple = TRUE,
                                   return_query_params = TRUE)
  expect_true(
    all(stringr::str_detect(setdiff(names(simp_query), names(simp)), "oncokb_query")))


  all_no_query <- annotate_mutations(oncokbR::blca_mutation[1:10, ],
                                 return_simple = FALSE,
                                 return_query_params = FALSE)

  all_with_query <- annotate_mutations(oncokbR::blca_mutation[1:10, ],
                                       return_simple = FALSE,
                                       return_query_params = TRUE)

  expect_true(
    all(stringr::str_detect(setdiff(names(all_with_query),
                                    names(all_no_query)), "oncokb_query")))

  expect_gt(length(names(simp_query)), length(names(simp)))
  expect_gt(length(names(all_with_query)), length(names(all_no_query)))

})
