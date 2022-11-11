

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

  annotate_tumor_type <- ("tumor_type" %in% names(cna))


  cna <- cna %>%
    mutate(hugo_symbol = as.character(.data$hugo_symbol)) %>%
    mutate(alteration = tolower(stringr::str_trim(as.character(.data$alteration))))

  levels_in_data <- names(table(cna$alteration))

  allowed_chr_levels <- c(
    "neutral" = "0",
    "DELETION" = "-2",
    "LOSS" = "-1.5",
    "LOSS" = "-1",
    "GAIN" = "1",
    "AMPLIFICATION" = "2"
  )

  all_allowed <- c(allowed_chr_levels, names(allowed_chr_levels))
  not_allowed <- levels_in_data[!levels_in_data %in% all_allowed]

  if(length(not_allowed) > 0) {
    cli::cli_abort(c("Unknown values in {.field alteration} field: {.val {not_allowed}}",
                     "Must be one of the following: {.val {all_allowed}}"))
  }


  suppressWarnings(
    cna <- cna %>%
      mutate(alteration_cleaned = forcats::fct_recode(.data$alteration, !!!allowed_chr_levels))
  )

  all_cna_oncokb <- cna %>%
    select(any_of(c("sample_id", "hugo_symbol", "alteration_cleaned", "tumor_type")))


  make_url <- function(sample_id, hugo_symbol, alteration_cleaned, tumor_type) {

    url <- glue::glue("https://www.oncokb.org/api/v1/annotate/copyNumberAlterations?hugoSymbol=",
                      "{hugo_symbol}&copyNameAlterationType={alteration_cleaned}&referenceGenome=GRCh37")

    if(annotate_tumor_type) {
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

  # Tumor Type - Remove Cols if None  ------------------------------------------
  if (!annotate_tumor_type) {
    all_cna_oncokb <- all_cna_oncokb %>%
      select(-contains("treatments"))
    cli::cli_alert_info("No {.val tumor_type} found in data. No treatment-level annotations will be returned.")
  }


  all_cna_oncokb

}

