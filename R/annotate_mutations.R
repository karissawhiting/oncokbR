

#' Annotate Mutations
#'
#' Annotate maf and maf-like files using information about protein change or HGVSg.
#' You can annotated with tumor types indicated (oncotree codes), or without.
#'
#' @param mutations a mutations file in MAF, or similar format
#' @param annotate_by Can indicate whether to annotate by "protein_change" or "hgvsg (see oncoKB API docs for more info).
#' "Default is `protein_change`
#'
#' @return an annotated mutations file
#' @export
#' @import dplyr
#'
#' @examples
#' mutations <- oncokbR::blca_mutation[1:10, ]
#'
#' # Annotate without tumor type -----
#' x <- annotate_mutations(mutations)
#'
#' # Annotate with tumor type -----

#'
#' # Annotate by HGVSg
#' z <- annotate_mutations(oncokbR::mutations_hgvsg, annotate_by =  "hgvsg")
#'
annotate_mutations <- function(mutations, annotate_by = c("protein_change", "hgvsg")) {

  mutations <- rename_columns(mutations)
  column_names <- colnames(mutations)
  annotate_by <- match.arg(annotate_by)

  # Protein Change -------------------------------------------------------------

  if(annotate_by == "protein_change") {

    # * Required Columns ------

    if(!("protein_pos_start" %in% column_names |
         "protein_pos_end" %in% column_names) &
       "protein_position" %in% column_names) {

        start_raw <- stringr::str_split_fixed(mutations$protein_position, "/", n=2)[,1]
        mutations$protein_pos_start <- stringr::str_split_fixed(start_raw, "-", n=2)[,1]
        mutations$protein_pos_end <- stringr::str_split_fixed(start_raw, "-", n=2)[,2]

        mutations <- mutations %>%
          mutate(protein_pos_end = case_when(
            protein_pos_end == "" ~ protein_pos_start,
          TRUE ~ protein_pos_end)) %>%
          mutate(across(c("protein_pos_start", "protein_pos_end"),
                        ~case_when(.x == "" ~ "NULL",
                                   TRUE ~ .x)))

      }

      required_cols_pro <- c("sample_id",
                             "hugo_symbol",
                             "hgv_sp_short",
                             "variant_classification",
                             "protein_pos_start",
                             "protein_pos_end")

      column_names <- colnames(mutations)
      which_missing <- required_cols_pro[which(!(required_cols_pro %in% column_names))]


      if(length(which_missing) > 0) {
        cli::cli_abort("The following required columns are missing in your mutations data: {.field {which_missing}}")
      }

      # * Check Variant Consequence  -----------

      variant_options <- unique(stats::na.omit(unlist(oncokbR::consequence_map)))
      variant_in_data <- unique(mutations$variant_classification)

      not_allowed <- stats::na.omit(variant_in_data[!(variant_in_data %in% variant_options)])

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

    }


  # HGVSG ----------------------------------------------------------------------

  # * Required Columns ------

  if(annotate_by == "hgvsg") {
    required_cols_pro <- c("sample_id", "hgv_sg")

    which_missing <- required_cols_pro[which(!(required_cols_pro %in% column_names))]


    if(length(which_missing) > 0) {
      cli::cli_abort("The following required columns are missing in your mutations data: {.field {which_missing}}")
    }
  }


  # Annotate Mutations  -------------------------------------------------------

  annotate_tumor_type <- ("tumor_type" %in% names(mutations))

  all_mut_oncokb <- switch(annotate_by,
                           "protein_change" =
                             .annotate_mutations_by_protein_change(mutations = mutations,
                                                                   annotate_tumor_type = annotate_tumor_type),
                           "hgvsg" = .annotate_mutations_by_hgvsg(mutations = mutations,
                                                                  annotate_tumor_type = annotate_tumor_type))

  # Clean Results  ----------------------------------------------------------

  all_mut_oncokb <- all_mut_oncokb %>%
    janitor::clean_names() %>%
    dplyr:: rename_with(~ stringr::str_remove(., "query_"),
                        .cols = starts_with("query_")) %>%
    select("sample_id", everything())

  # If no tumor type, remove those columns
  if (!annotate_tumor_type) {
    all_mut_oncokb <- all_mut_oncokb %>%
      select(-contains("treatments"))
    cli::cli_alert_info("No {.val tumor_type} found in data. No treatment-level annotations will be returned.")
  }

  all_mut_oncokb

}

# Sub Annotators By Data Type ----------------------------------------------------

#' Annotate Mutations By Protein Change
#'
#' @param mutations
#'
#' @return a dataframe of annotated mutations
#' @export
#'
#' @keywords internal
#' @examples
#' mutations <- rename_columns(oncokbR::blca_mutation[1:10, ]) %>%
#'     dplyr::mutate(consequence_final_coding = "any")
#' .annotate_mutations_by_protein_change(mutations, annotate_tumor_type = FALSE)
#'
.annotate_mutations_by_protein_change <- function(mutations, annotate_tumor_type) {

  all_mut_oncokb <- mutations %>%
    select(any_of(c("sample_id", "hugo_symbol", "hgv_sp_short", "consequence_final_coding",
                    "protein_pos_start", "protein_pos_end", "tumor_type")))

  make_url <- function(sample_id, hugo_symbol, hgv_sp_short,
                       consequence_final_coding,
                       protein_pos_start, protein_pos_end, tumor_type) {

    url <- glue::glue("https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?hugoSymbol=",
                      "{hugo_symbol}&alteration={hgv_sp_short}&referenceGenome=GRCh37&consequence=",
                      "{consequence_final_coding}&proteinStart={protein_pos_start}&proteinEnd={protein_pos_end}")


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
    # parsed <- parsed %>% discard("query")
    parsed <-unlist(parsed, recursive=TRUE) %>%
      tibble::enframe() %>%
      tidyr::pivot_wider(names_from = .data$name,
                         values_fn = function(x) paste(x, collapse=","))

    parsed$sample_id <- sample_id
    parsed

  }

  purrr::pmap_df(all_mut_oncokb,  make_url)


 }



#' Annotate Mutations By HGVSg
#'
#' @param mutations
#'
#' @return a dataframe of annotated mutations
#' @export
#'
#' @keywords internal
#' @examples
#'
#'
#' mutations <- rename_columns(oncokbR::mutations_hgvsg)
#' .annotate_mutations_by_hgvsg(mutations, annotate_tumor_type = FALSE)
#'
.annotate_mutations_by_hgvsg <- function(mutations, annotate_tumor_type) {

  all_mut_oncokb <- mutations %>%
    select(any_of(c("sample_id", "hgv_sg", "tumor_type")))

  make_url <- function(sample_id, hgv_sg, tumor_type) {

    url <- glue::glue("https://www.oncokb.org/api/v1/annotate/mutations/byHGVSg?hgvsg=",
                      "{hgv_sg}&referenceGenome=GRCh37")

    if(annotate_tumor_type) {
      url <- glue::glue(url, "&tumorType={tumor_type}")
    }

    url <- utils::URLencode(url)

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

  purrr::pmap_df(all_mut_oncokb,  make_url)


}



