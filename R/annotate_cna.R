

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
                         annotate_tumor_type = NULL,
                         return_simple = TRUE,
                         return_query_params = FALSE,
                         reference_genome = "GRCh37",
                         token = get_oncokb_token()) {

  # Check Data --------------------------------------------------------------

  # standardize column names
  cna <- rename_columns(cna)

  # check for tumor type
  annotate_tumor_type <- ("tumor_type" %in% names(cna))

  # Check required columns & data types ---------------------------------------
  .check_required_cols(data = cna,
                       required_cols = c("sample_id", "hugo_symbol", "alteration"),
                       data_name = "cna")

  # Data Pre-processing -----------------------------------------------------

  cna <- cna %>%
    mutate(hugo_symbol = as.character(.data$hugo_symbol)) %>%
    mutate(alteration = tolower(stringr::str_trim(as.character(.data$alteration))))

  # * Recode Alterations ------

  levels_in_data <- names(table(cna$alteration))

  cna <- cna %>%
    mutate(alteration = recode_cna(.data$alteration)) %>%
    mutate(alteration = toupper(.data$alteration))


  # Annotate CNA ------------------------------------------------------------

  # * Check Token ------------
  validate_oncokb_token(token = token)

  # * Check Parameters ------------

  # `reference_genome`
  if (!inherits(referenceGenome, "character")) {
    stop("`referenceGenome` must be a character")
  }

  # * Make request -------

  # create event index
  cna <- mutate(cna, event_index = 1:nrow(cna))

  cna_select <- cna %>%
    select(any_of(c("event_index",
                    "sample_id",
                    "hugo_symbol",
                    "alteration",
                    "tumor_type")))

  # If no tumor type, give NAs for call
  if(!annotate_tumor_type) {
    cna_select <- cna_select %>%
      mutate(tumor_type = NA_character_)
  }

  # Todo: This could be cleaned up or even separated into own function
  requests <- pmap(cna_select, ~oncokb_api_request(
    event_index = ..1,
    sample_id = ..2,
    resource = "annotate/copyNumberAlterations",
    hugoSymbol = ..3,
    copyNameAlterationType = ..4,
    tumor_type = ..5,
    referenceGenome = reference_genome,
    token = token
  ))

  list_of_responses <- map(requests, function(request) {
    req_perform(request) |>
      httr2::resp_body_string()})

  all_cna_oncokb <- map_df(list_of_responses, parse_responses)
  all_cna_oncokb <- bind_cols(cna_select[, c('event_index', 'sample_id')],
                              all_cna_oncokb)

  # Clean Results  ----------------------------------------------------------

  all_cna_oncokb <- .clean_query_results(
    query_result = all_cna_oncokb,
    return_simple = return_simple,
    return_query_params = return_query_params,
    original_data = cna)

  # * Tumor Type - Remove Col if None  ------------------------------------------

  all_cna_oncokb <- .tumor_type_warning(
    annotate_tumor_type = annotate_tumor_type,
    data = all_cna_oncokb)

  return(all_cna_oncokb)

}

