

test_that("Message thrown when no tumor_type", {

  # mutations -----
  expect_message(mut <-
    annotate_mutations(oncokbR::blca_mutation[1:10, ]))

  expect_no_message(mut_tum <-
    annotate_mutations(
      mutate(oncokbR::blca_mutation[1:10, ], tumor_type = "BLCA")))

  # cna -------
  expect_message(cna <- annotate_cna(oncokbR::blca_cna[1:10, ]))

  expect_no_message(
    cna_tum <- annotate_cna(
      mutate(oncokbR::blca_cna[1:10, ], tumor_type = "BLCA")))

  # sv -------
  expect_message(
    sv <- annotate_sv(oncokbR::blca_sv[1:10, ]))

  expect_no_message(
    sv_tum <- annotate_sv(
      mutate(oncokbR::blca_sv[1:10, ], tumor_type = "BLCA")))

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
  expect_equal(length(col_diff), 47)
  expect_equal(sum(stringr::str_detect(col_diff, "treatments")),
               46)

  # cna -------
  cna_tum <- annotate_cna(
      mutate(oncokbR::blca_cna[1:10, ], tumor_type = "BLCA"),
      return_simple = FALSE)

  cna <- annotate_cna(
    oncokbR::blca_cna[1:10, ],
    return_simple = FALSE)

  col_diff <- setdiff(names(cna_tum), names(cna))
  expect_equal(length(col_diff), 70)
  expect_equal(sum(stringr::str_detect(col_diff, "treatments")),
               69)

  # sv -------
  sv_tum <- annotate_sv(
      mutate(oncokbR::blca_sv[1:10, ], tumor_type = "BLCA"),
      return_simple = FALSE)

  sv <- annotate_sv(
    oncokbR::blca_sv[1:10, ],
    return_simple = FALSE)

  col_diff <- setdiff(names(sv_tum), names(sv))
  expect_equal(length(col_diff), 24)
  expect_equal(sum(stringr::str_detect(col_diff, "treatments")),
               22)
})
