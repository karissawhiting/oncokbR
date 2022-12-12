

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

  # Check required columns & data types ---------------------------------------
  required_cols <- c("sample_id", "hugo_symbol", "alteration")
  column_names <- colnames(cna)

  which_missing <- required_cols[which(!(required_cols %in% column_names))]

  if(length(which_missing) > 0) {
    cli::cli_abort("The following required columns are missing in your mutations data: {.field {which_missing}}")
  }

  # Clean Data --------------------------------------------------------------
  cna <- cna %>%
    mutate(hugo_symbol = as.character(.data$hugo_symbol)) %>%
    mutate(alteration = tolower(stringr::str_trim(as.character(.data$alteration))))

  levels_in_data <- names(table(cna$alteration))

  # recode alterations
  cna <- cna %>%
    mutate(alteration = recode_cna(.data$alteration)) %>%
    mutate(alteration = toupper(.data$alteration))

  all_cna_oncokb <- cna %>%
    select(any_of(c("sample_id", "hugo_symbol", "alteration", "tumor_type")))


  make_url <- function(sample_id, hugo_symbol, alteration, tumor_type) {

    url <- glue::glue("https://www.oncokb.org/api/v1/annotate/copyNumberAlterations?hugoSymbol=",
                      "{hugo_symbol}&copyNameAlterationType={alteration}&referenceGenome=GRCh37")

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

  # Clean Results  ----------------------------------------------------------

  all_cna_oncokb <- all_cna_oncokb %>%
    janitor::clean_names() %>%
    dplyr:: rename_with(~ stringr::str_remove(., "query_"),
                        .cols = starts_with("query_")) %>%
    select("sample_id", everything())

  # Tumor Type - Remove Cols if None  ------------------------------------------
  if (!annotate_tumor_type) {
    all_cna_oncokb <- all_cna_oncokb %>%
      select(-contains("treatments"))
    cli::cli_alert_info("No {.val tumor_type} found in data. No treatment-level annotations will be returned.")
  }


  all_cna_oncokb

}

