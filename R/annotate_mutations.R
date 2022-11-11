

#' Annotate Mutations
#'
#' Annotate maf and maf-like files using information about protein change or HGVSg.
#' You can annotated with tumor types indicated (oncotree codes), or without.
#'
#' @param mutations a mutations file in MAF, or similar format
#' @param annotate_by Can indicate whether to annotate by "protein_change" or "hgvsg (see oncoKB API docs for more info).
#' If none specified, a guess will be made according to columns in your data.
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
#' mutations$tumor_type = "BLCA"
#' y <- annotate_mutations(mutations)
#'
#' # Annotate by HGVSg
#' z <- annotate_mutations(oncokbR::mutations_hgvsg)
#'
annotate_mutations <- function(mutations, annotate_by = NULL) {

  ## ADD IN WHICH COLUMNS ARE REQUIRED

  mutations <- rename_columns(mutations)

  annotate_by <- switch(!missing(annotate_by),
         match.arg(annotate_by, c("protein_change", "hgvsg")))

  annotate_by <- annotate_by %||% {
    case_when(
      all(c("hgv_sp_short", "protein_pos_start", "protein_pos_end") %in% names(mutations)) ~ "protein_change",
    "hgv_sg" %in% names(mutations) ~ "hgvsg")
    }


  annotate_tumor_type <- ("tumor_type" %in% names(mutations))

  # Check Variant Conseqeunce  -------------------------------------------------
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


  # Annotate Mutations  -------------------------------------------------------
  all_mut_oncokb <- switch(annotate_by,
                           "protein_change" =
                             .annotate_mutations_by_protein_change(mutations = mutations,
                                                                   annotate_tumor_type = annotate_tumor_type),
                           "hgvsg" = .annotate_mutations_by_hgvsg(mutations = mutations,
                                                                  annotate_tumor_type = annotate_tumor_type))

  all_mut_oncokb <- all_mut_oncokb %>%
    janitor::clean_names() %>%
    dplyr:: rename_with(~ stringr::str_remove(., "query_"), .cols = starts_with("query_")) %>%
    select("sample_id", everything())

  # Tumor Type - Remove Cols if None  ------------------------------------------
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
#' mutations <- rename_columns(oncokbR::blca_mutation[1:10, ])
#' .annotate_mutations_by_protein_change(mutations)
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
#' .annotate_mutations_by_hgvsg(mutations)
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

    url <- URLencode(url)

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



