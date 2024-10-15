

test_that("Message thrown when no tumor_type", {

  # mutations -----
  expect_message(mut2 <-
    annotate_mutations(oncokbR::blca_mutation[1:10, ]),
    regexp = "No treatment-level annotations will be returned")

  expect_no_message(mut_tum <-
    annotate_mutations(
      mutate(oncokbR::blca_mutation[1:10, ], tumor_type = "BLCA")),
    message = "No treatment-level annotations will be returned")

  # cna -------
  expect_message(cna <- annotate_cna(oncokbR::blca_cna[1:10, ], return_simple = FALSE, return_query_params = TRUE),
                 regexp = "No treatment-level annotations will be returned")

  expect_no_message(
    cna_tum <- annotate_cna(
      mutate(oncokbR::blca_cna[1:10, ], tumor_type = "BLCA"), return_simple = FALSE, return_query_params = TRUE),
    message = "No treatment-level annotations will be returned")

  # sv -------
  expect_message(
    sv <- annotate_sv(oncokbR::blca_sv[1:10, ]),
    regexp = "No treatment-level annotations will be returned")

  expect_no_message(
    sv_tum <- annotate_sv(
      mutate(oncokbR::blca_sv[1:10, ], tumor_type = "BLCA"
             )),
    message = "No treatment-level annotations will be returned")

  # Data different with tumor type ---
  expect_true(na.omit(mut$oncokb_highest_sensitive_level) !=
                 na.omit(mut_tum$oncokb_highest_sensitive_level))

  expect_true(cna$oncokb_highest_sensitive_level[6] !=
                cna_tum$oncokb_highest_sensitive_level[6])

  expect_true(na.omit(sv$oncokb_highest_sensitive_level) !=
                na.omit(sv_tum$oncokb_highest_sensitive_level))

})


test_that("Correct columns returned with tumor_type", {

  # mutations -----
  mut_tum <- annotate_mutations(
      mutate(oncokbR::blca_mutation[1:10, ], tumor_type = "BLCA"),
      return_simple = FALSE)

  mut <- annotate_mutations(
    oncokbR::blca_mutation[1:10, ],
    return_simple = FALSE)

  col_diff <- setdiff(names(mut_tum), names(mut))
  expect_true(
    all(stringr::str_detect(col_diff[2:length(col_diff)], 'treatments')))

  # cna -------
  cna_tum <- annotate_cna(
      mutate(oncokbR::blca_cna[1:10, ], tumor_type = "BLCA"),
      return_simple = FALSE)

  cna <- annotate_cna(
    oncokbR::blca_cna[1:10, ],
    return_simple = FALSE)

  col_diff <- setdiff(names(cna_tum), names(cna))
  expect_true(
    all(stringr::str_detect(col_diff[2:length(col_diff)], 'treatments')))


  # sv -------
  sv_tum <- annotate_sv(
      mutate(oncokbR::blca_sv[1:10, ], tumor_type = "BLCA"),
      return_simple = FALSE)

  sv <- annotate_sv(
    oncokbR::blca_sv[1:10, ],
    return_simple = FALSE)

  col_diff <- setdiff(names(sv_tum), names(sv))
  col_diff <- col_diff[!stringr::str_detect(col_diff, "tumor_type")]
  expect_true(
    all(stringr::str_detect(col_diff[2:length(col_diff)], 'treatments')))

})
