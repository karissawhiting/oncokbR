

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
#'
annotate_sv <- function(sv,
                        return_simple = TRUE,
                        return_query_params = FALSE) {

  sv <- rename_columns(sv)
  annotate_tumor_type <- ("tumor_type" %in% names(sv))

  # Check required columns & data types ---------------------------------------
  required_cols_sv <- c("sample_id", "site_1_hugo_symbol", "site_2_hugo_symbol",
                        "variant_class")

  .check_required_cols(data = sv,
                       required_cols = required_cols_sv,
                       data_name = "sv")

  # Clean Data --------------------------------------------------------------
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

  # Clean Variant Class -----------------------------------------------------

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


  # Annotate SV ----------------------------------------------------------------

  sv <- mutate(sv, index = 1:nrow(sv))

  all_sv_oncokb_raw <- sv %>%
    select(any_of(c("index", "sample_id",
                    "site_1_hugo_symbol", "site_2_hugo_symbol", "structural_variant_type", "is_functional", "tumor_type")))

  make_url <- function(index, sample_id, site_1_hugo_symbol, site_2_hugo_symbol,
                       structural_variant_type, is_functional, tumor_type) {

    url <- glue::glue("https://www.oncokb.org/api/v1/annotate/structuralVariants?hugoSymbolA=",
                      "{site_1_hugo_symbol}&hugoSymbolB={site_2_hugo_symbol}&structuralVariantType=",
                      "{structural_variant_type}&isFunctionalFusion={is_functional}&referenceGenome=GRCh37")


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
      tidyr::pivot_wider(names_from = "name",
                         values_fn = function(x) paste(x, collapse=","))

    parsed$sample_id <- sample_id
    parsed$index <- index

    parsed
  }


  all_sv_oncokb <- purrr::pmap_df(all_sv_oncokb_raw,  make_url)

  # Clean Results  ----------------------------------------------------------

  all_sv_oncokb <- all_sv_oncokb %>%
    janitor::clean_names() %>%
    dplyr:: rename_with(~ stringr::str_remove(., "query_"), .cols = starts_with("query_")) %>%
    rename("query_hugo_symbol" = "hugo_symbol") %>%
    select("sample_id", everything())

  all_sv_oncokb <- .clean_query_results(
    query_result = all_sv_oncokb,
    return_simple = return_simple,
    return_query_params = return_query_params,
    original_data = sv)

  # Tumor Type - Remove Cols if None  ------------------------------------------

  all_sv_oncokb <- .tumor_type_warning(
    annotate_tumor_type = annotate_tumor_type,
    data = all_sv_oncokb)

  return(all_sv_oncokb)


}

