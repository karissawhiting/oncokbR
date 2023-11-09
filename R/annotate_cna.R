

#' Annotate CNA
#'
#' @param cna a cna file in long format (similar to maf file)
#' @inheritParams annotate_mutations
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
annotate_cna <- function(cna,
                         return_simple = TRUE,
                         return_query_params = FALSE) {

  cna <- rename_columns(cna)
  annotate_tumor_type <- ("tumor_type" %in% names(cna))

  # Check required columns & data types ---------------------------------------
  required_cols_cna <- c("sample_id", "hugo_symbol", "alteration")

  .check_required_cols(data = cna,
                       required_cols = required_cols_cna,
                       data_name = "cna")

  # Clean Data --------------------------------------------------------------
  cna <- cna %>%
    mutate(hugo_symbol = as.character(.data$hugo_symbol)) %>%
    mutate(alteration = tolower(stringr::str_trim(as.character(.data$alteration))))

  levels_in_data <- names(table(cna$alteration))

  # recode alterations
  cna <- cna %>%
    mutate(alteration = recode_cna(.data$alteration)) %>%
    mutate(alteration = toupper(.data$alteration))



# Annotate CNA ------------------------------------------------------------

  cna <- mutate(cna, index = 1:nrow(cna))

  all_cna_oncokb <- cna %>%
    select(any_of(c("index", "sample_id", "hugo_symbol", "alteration", "tumor_type")))


  make_url <- function(index, sample_id, hugo_symbol,
                       alteration, tumor_type) {

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

    parsed <- unlist(parsed, recursive=TRUE) %>%
      tibble::enframe() %>%
      tidyr::pivot_wider(names_from = "name",
                  values_fn = function(x) paste(x, collapse=","))

    parsed$sample_id <- sample_id
    parsed$index <- index
    parsed
  }

  all_cna_oncokb <- purrr::pmap_df(all_cna_oncokb,  make_url)

  # Clean Results  ----------------------------------------------------------

  all_cna_oncokb <- .clean_query_results(
    query_result = all_cna_oncokb,
    return_simple = return_simple,
    return_query_params = return_query_params,
    original_data = cna)

  return(all_cna_oncokb)

}

