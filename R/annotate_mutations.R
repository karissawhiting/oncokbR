

#' Annotate Mutations
#'
#' @param mutations a mutations file in MAF, or similar format
#'
#' @return an annotated mutations file
#' @export
#' @import dplyr
#'
#' @examples
#' ex_mut <- blca_mutation[1:50, ]
#'
#' # Annotate without tumor type -----
#' x <- annotate_mutations(ex_mut)
#'
#' # Annotate with tumor type -----
#' ex_mut$tumor_type = "BLCA"
#' y <- annotate_mutations(ex_mut)
#'
annotate_mutations <- function(mutations) {

  mutations <- rename_columns(mutations)

  variant_options <- unique(stats::na.omit(unlist(oncokbR::consequence_map)))
  variant_in_data <- unique(mutations$variant_classification)

  not_allowed <- variant_in_data[!(variant_in_data %in% variant_options)]

  # Maybe turn into warning
  if(length(not_allowed) > 0) {
    cli::cli_abort("The following variant classification levels are not recognized: {.code {not_allowed}}.
                   Please remove or recode these to continue (see {.code oncokbR::consequence_map} for allowed values)")
  }

  consequence_vec <- tibble::deframe(select(
    oncokbR::consequence_map, "consequence_final_coding", "variant_classification"))

  suppressWarnings(
    mutations <- mutations %>%
      mutate(consequence_final_coding = forcats::fct_recode(.data$variant_classification, !!!consequence_vec))
  )

  all_mut_oncokb <- mutations %>%
    select(any_of(c("sample_id", "hugo_symbol", "hgv_sp_short", "consequence_final_coding",
           "protein_pos_start", "protein_pos_end", "tumor_type")))

  make_url <- function(sample_id, hugo_symbol, hgv_sp_short,
                       consequence_final_coding,
                       protein_pos_start, protein_pos_end, tumor_type) {

    url <- glue::glue("https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?hugoSymbol=",
                      "{hugo_symbol}&alteration={hgv_sp_short}&referenceGenome=GRCh37&consequence=",
                      "{consequence_final_coding}&proteinStart={protein_pos_start}&proteinEnd={protein_pos_end}")


    if(("tumor_type" %in% names(mutations))) {
      url <- glue::glue(url, "&tumorType={tumor_type}")
    }

    token <- Sys.getenv('ONCOKB_TOKEN')
    resp <- httr::GET(url,
                      httr::add_headers(Authorization = paste("Bearer ",
                                                              token,
                                                              sep = "")))

    parsed <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                                 flatten = TRUE,
                                 simplifyVector = TRUE)
    # parsed <- parsed %>% discard("query")
    parsed <-unlist(parsed, recursive=TRUE) %>%
      tibble::enframe() %>%
      tidyr::pivot_wider(names_from = .data$name,
                  values_fn = function(x) paste(x, collapse=","))

    parsed$sample_id <- sample_id
    parsed

  }


  all_mut_oncokb <- purrr::pmap_df(all_mut_oncokb,  make_url)

  all_mut_oncokb <- all_mut_oncokb %>%
    janitor::clean_names() %>%
    dplyr:: rename_with(~ stringr::str_remove(., "query_"), .cols = starts_with("query_")) %>%
    select("sample_id", everything())

  all_mut_oncokb

}

