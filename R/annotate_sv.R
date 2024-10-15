

#' Annotate Structural Variants
#'
#' @param sv a sv file in long format (similar to maf file)
#' @inheritParams annotate_mutations
#' @return an annotated sv file
#' @export
#' @import dplyr
#'
#' @examples
#' sv <- blca_sv[1:10,]
#'
#' x <- annotate_sv(sv = sv)
#' x <- annotate_sv(sv = mutate(oncokbR::blca_sv[1:10, ], tumor_type = "BLCA"))
annotate_sv <- function(sv,
                        annotate_tumor_type = NULL,
                        return_simple = TRUE,
                        return_query_params = FALSE,
                        reference_genome = "GRCh37",
                        token = get_oncokb_token()) {

  # Check Data --------------------------------------------------------------

  # standardize column names
  sv <- rename_columns(sv)

  # check for tumor type
  annotate_tumor_type <- ("tumor_type" %in% names(sv))

  # check other columns
  .check_required_cols(data = sv,
                       required_cols = c("sample_id",
                                         "site_1_hugo_symbol",
                                         "site_2_hugo_symbol",
                                         "variant_class"),
                       data_name = "sv")

  # Data Pre-processing -----------------------------------------------------

  sv <- sv %>%
    mutate(site_1_hugo_symbol = as.character(.data$site_1_hugo_symbol)) %>%
    mutate(site_2_hugo_symbol = as.character(.data$site_2_hugo_symbol)) %>%
    mutate(structural_variant_type =
             stringr::str_trim(as.character(.data$variant_class))) %>%
    mutate(structural_variant_type =
             case_when(
               .data$structural_variant_type %in% c("NA") |
                 .data$structural_variant_type %in% c("") |
                 is.na(.data$structural_variant_type) ~ "UNKNOWN",
           TRUE ~ .data$structural_variant_type)) %>%
    mutate(structural_variant_type = toupper(.data$structural_variant_type))


  # * Determine Functional Variants ---------------
  # Assume all structural variants are functional (this mirrors behavior in
  # https://github.com/oncokb/oncokb-annotator/blob/47e4a158ee843ead75445982532eb149db7f3106/AnnotatorCore.py#L1506)

  # intragenic not counted as functional (as per python annotator)
  if(!("is_functional" %in% names(sv))) {
    sv <- sv %>%
      mutate(is_functional =
               case_when(
                 site_1_hugo_symbol == site_1_hugo_symbol ~ "false",
                 TRUE ~ "true"))
    }

  # * Clean Variant Class ---------------

  levels_in_data <- names(table(sv$structural_variant_type))

  allowed_chr_levels <- c("DELETION",
                          "TRANSLOCATION",
                          "DUPLICATION",
                          "INSERTION",
                          "INVERSION",
                          "FUSION",
                          "UNKNOWN")

  all_allowed <- c(allowed_chr_levels, names(allowed_chr_levels))
  not_allowed <- levels_in_data[!levels_in_data %in% all_allowed]

  if(length(not_allowed) > 0) {
    cli::cli_abort(c("Unknown values in {.field variant_class} field: {.val {not_allowed}}",
                     "Must be one of the following: {.val {all_allowed}}"))
  }



  # Annotate SV  ---------------------------------------------------------

  # * Check Token ------------
  validate_oncokb_token(token = token)

  # * Check Parameters ------------

  # `reference_genome``
  if (!inherits(reference_genome, "character")) {
    stop("`reference_genome` must be a character")
  }


  # * Make request -------

  # create index
  sv <- sv %>%
    mutate(event_index = 1:nrow(.))

  sv_select <- sv %>%
    dplyr::select(any_of(c(
      "event_index",
      "sample_id",
      "site_1_hugo_symbol",
      "site_2_hugo_symbol",
      "structural_variant_type",
      "is_functional",
      "tumor_type")))

  # If no tumor type, give NAs for call
  if(!annotate_tumor_type) {
    sv_select <- sv_select %>%
      mutate(tumor_type = NA_character_)
  }

  # Todo: This could be cleaned up or even separated into own function
  requests <- pmap(sv_select, ~oncokb_api_request(
    event_index = ..1,
    sample_id = ..2,
    resource = "annotate/structuralVariants",
    hugoSymbolA = ..3,
    hugoSymbolB = ..4,
    structuralVariantType = ..5,
    isFunctionalFusion = ..6,
    tumorType = ..7,
    token = token,
    referenceGenome = reference_genome
  ))


  list_of_responses <- map(requests, function(request) {
    req_perform(request) |>
      httr2::resp_body_string()})

  all_sv_oncokb <- map_df(list_of_responses, parse_responses)
  all_sv_oncokb <- bind_cols(sv_select[, c('event_index', 'sample_id')],
                             all_sv_oncokb)


  # Clean Results  ----------------------------------------------------------

  all_sv_oncokb <- all_sv_oncokb %>%
    janitor::clean_names() %>%
    dplyr:: rename_with(~ stringr::str_remove(., "query_"), .cols = starts_with("query_")) %>%
    rename("query_hugo_symbol" = "hugo_symbol") %>%
    select("event_index", "sample_id", everything())

  all_sv_oncokb <- .clean_query_results(
    query_result = all_sv_oncokb,
    return_simple = return_simple,
    return_query_params = return_query_params,
    original_data = sv)

  # * Tumor Type - Remove col if None  -------

  all_sv_oncokb <- .tumor_type_warning(
    annotate_tumor_type = annotate_tumor_type,
    data = all_sv_oncokb)

  return(all_sv_oncokb)


}
