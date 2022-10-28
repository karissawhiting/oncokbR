

#' Annotate CNA
#'
#' @param cna a cna file in long format (similar to maf file)
#'
#' @return an annotated cna file
#' @export
#' @import dplyr
#'
#' @examples
#' ex_cna <- blca_cna[1:10,]
#' ex_cna$tumor_type = "BLCA"
#'
#' x <- annotate_cna(ex_cna)
#'
annotate_cna <- function(cna) {

  cna <- rename_columns(cna)

  if(!("tumor_type" %in% names(cna))) {
    cli::cli_abort("Your data must have a column with {.val tumor_type}")
  }


  cna <- switch(!is.null(cna), .recode_cna_alterations(cna))

  cna <- cna %>%
    mutate(alteration_cleaned = alteration)

  all_cna_oncokb <- cna %>%
    select("sample_id", "hugo_symbol", "alteration_cleaned", "tumor_type")

  make_url <- function(sample_id, hugo_symbol, alteration_cleaned, tumor_type) {

    url <- glue::glue("https://www.oncokb.org/api/v1/annotate/copyNumberAlterations?hugoSymbol=",
                      "{hugo_symbol}&copyNameAlterationType={alteration_cleaned}&referenceGenome=GRCh37&",
                      "tumorType={tumor_type}")

    token <- Sys.getenv('ONCOKB_TOKEN')
    resp <- httr::GET(url,
                      httr::add_headers(Authorization = paste("Bearer ",
                                                              token,
                                                              sep = "")))

    parsed <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                                 flatten = TRUE,
                                 simplifyVector = TRUE)

    parsed <-unlist(parsed, recursive=TRUE) %>%
      tibble::enframe() %>%
      tidyr::pivot_wider(names_from = .data$name,
                  values_fn = function(x) paste(x, collapse=","))

    parsed$sample_id <- sample_id
    parsed
  }


  all_cna_oncokb <- purrr::pmap_df(all_cna_oncokb,  make_url)

  all_cna_oncokb <- all_cna_oncokb %>%
    janitor::clean_names() %>%
    select(-contains("query_")) %>%
    select("sample_id", everything())

  all_cna_oncokb

}

