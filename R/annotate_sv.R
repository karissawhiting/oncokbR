

#' Annotate CNA
#'
#' @param cna a cna file in long format (similar to maf file)
#'
#' @return an annotated cna file
#' @export
#' @import dplyr
#'
#' @examples
#' ex_sv <- blca_sv[1:10,]
#'
#' x <- annotate_sv(sv = ex_sv)
#'
annotate_sv <- function(sv) {

  sv <- rename_columns(sv)

  # if(!("tumor_type" %in% names(cna))) {
  #   cli::cli_abort("Your data must have a column with {.val tumor_type}")
  # }

  sv <- sv %>%
    mutate(hugo_symbol = as.character(.data$site_1_hugo_symbol)) %>%
    mutate(hugo_symbol = as.character(.data$site_2_hugo_symbol)) %>%
    mutate(structural_variant_type =
             toupper(stringr::str_trim(as.character(.data$variant_class)))) %>%
    mutate(structural_variant_type = case_when(
      (is.na(.data$variant_class) | structural_variant_type == "NA") ~ "UNKNOWN",
      TRUE ~ structural_variant_type
    ))


  # Assume all structural variants are functions (this mirrors behavior in
  # https://github.com/oncokb/oncokb-annotator/blob/47e4a158ee843ead75445982532eb149db7f3106/AnnotatorCore.py#L1506)
  if(!("is_functional" %in% names(sv))) {
    sv <- sv %>%
      mutate(is_functional = "true")

  }

  levels_in_data <- names(table(sv$variant_class))

  # explicitely code NAs
  if("NA" %in% levels_in_data) {
    sv <- sv %>%
      mutate(variant_class =
        case_when(
          (.data$variant_class %in% "NA" |
            is.na(.data$variant_class)) ~ "UNKNOWN",
          TRUE ~ .data$variant_class))
  }

  levels_in_data <- names(table(sv$variant_class))

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
    cli::cli_abort(c("Unknown values in {.field alteration} field: {.val {not_allowed}}",
                     "Must be one of the following: {.val {all_allowed}}"))
  }


  all_sv_oncokb <- sv %>%
    select("sample_id", "site_1_hugo_symbol", "site_2_hugo_symbol", "structural_variant_type", "is_functional")

  make_url <- function(sample_id, site_1_hugo_symbol, site_2_hugo_symbol,
                       structural_variant_type, is_functional) {

    url <- glue::glue("https://www.oncokb.org/api/v1/annotate/structuralVariants?hugoSymbolA=",
                      "{site_1_hugo_symbol}&hugoSymbolB={site_2_hugo_symbol}&structuralVariantType=",
                      "{structural_variant_type}&isFunctionalFusion={is_functional}&referenceGenome=GRCh37")


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


  all_sv_oncokb <- purrr::pmap_df(all_sv_oncokb,  make_url)

  all_sv_oncokb <- all_sv_oncokb %>%
    janitor::clean_names() %>%
    select(-contains("query_")) %>%
    select("sample_id", everything())

  all_sv_oncokb

}
