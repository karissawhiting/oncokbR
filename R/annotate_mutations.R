#' Annotate Mutations
#'
#' Annotate maf and maf-like files using information about protein change or HGVSg.
#' You can annotated with tumor types indicated (oncotree codes), or without.
#'
#' @param mutations a mutations file in MAF, or similar format
#' @param annotate_by Can indicate whether to annotate by "protein_change" or "hgvsg (see oncoKB API docs for more info).
#' "Default is `protein_change`
#' @param return_simple Default is `TRUE` where only a set of the most common columns are returned from oncoKB annotator
#' see `oncokbR::output_dictionary` for more information on what is returned. If `FALSE` all raw columns are returned from API.
#' @param return_query_params If `TRUE`, the specific parameters used to query the API are returned in new columns.
#' This can be useful for troubleshooting the annotator. Default is `FALSE`.
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
#' z <- annotate_mutations(mutations = oncokbR::mutations_hgvsg, annotate_by =  "hgvsg")
#' z <- annotate_mutations(mutations = oncokbR::blca_mutation, annotate_by =  "protein_change")
#'

# Prep and Check Data --------------------------------------------------------
annotate_mutations <- function(mutations,
                               annotate_by = c("protein_change", "hgvsg"),
                               annotate_tumor_type = NULL,
                               return_simple = TRUE,
                               return_query_params = FALSE,
                               reference_genome = "GRCh37",
                               token = get_oncokb_token()) {
  # Check Annotation Method ---------------------------------------------------
  annotate_by <- match.arg(annotate_by)

  # Check Data --------------------------------------------------------------

  # standardize column names
  mutations <- rename_columns(mutations)
  column_names <- colnames(mutations)

  # check for tumor type
  annotate_tumor_type <- ("tumor_type" %in% column_names)

  # check other columns
  required_cols <- switch(annotate_by,
    "hgvsg" = c("sample_id", "hgv_sg"),
    "protein_change" = c(
      "sample_id",
      "hugo_symbol",
      "hgv_sp_short",
      "variant_classification",
      "protein_pos_start",
      "protein_pos_end"
    )
  )


  .check_required_cols(
    data = mutations,
    required_cols = required_cols,
    data_name = "mutations"
  )


  # Check Parameters ----------------------------------------------------------

  # * Reference Genome
  if (!inherits(reference_genome, "character")) {
    stop("`reference_genome` must be a character")
  }

  # Annotate Mutations --------------------------------------------------------

  # Set Resource
  resource <- switch(annotate_by,
    "hgvsg" = "annotate/mutations/byHGVSg",
    "protein_change" = "annotate/mutations/byProteinChange"
  )

  # Check Token
  validate_oncokb_token(token = token)

  # create event index
  mutations <- mutate(mutations, event_index = 1:nrow(mutations))

  # subset data to only necessary columns
  mutations_select <- mutations %>%
    select(any_of(c("event_index", required_cols, "tumor_type")))

  # * If HGVSG ------------------------------------------------

  if (annotate_by == "hgvsg") {
    if (!inherits(mutations_select$hgv_sg, "character")) {
      stop("`hgv_sg` must be a character")
    }

    requests <- map(mutations_select$hgv_sg, ~ oncokb_api_request(
      resource = resource,
      hgvsg = .x,
      referenceGenome = reference_genome,
      token = token
    ))

    list_of_responses <- map(requests, function(request) {
      req_perform(request) |>
        httr2::resp_body_string()
    })

    all_mutations_oncokb <- map_df(list_of_responses, parse_responses)

    all_mutations_oncokb <- bind_cols(
      mutations_select[, c("event_index", "sample_id")],
      all_mutations_oncokb
    )

    #**  Check errors ---------
    if ("message" %in% colnames(all_mutations_oncokb)) {
      messages <- stats::na.omit(all_mutations_oncokb$message)

      if (length(messages) > 0) {
        cli::cli_warn("{.val {length(messages)}} obervation(s) came back with the following information: {.code {unique(messages)}}
                        \nThese may be due to missing/NA data in columns used for annotation. To further check these values, see the {.col message} and {.col status}
                      columns of your resulting dataframe (e.g. {.code results[which(!is.na(results$message)), c('message', 'status')] })")
      }
    }
  }

  # * If Protein Change ------------------------------------------------

  if (annotate_by == "protein_change") {
    # create a protein start/end column from protein position if needed
    mutations_select <- .check_protein_start_end(column_names, mutations_select)


    # ** Check Variant Consequence  -----------
    .check_consequence(mutations_select$variant_classification)

    # recode vector as needed (consider combining with warning)
    consequence_vec <- tibble::deframe(select(
      oncokbR::consequence_map, "consequence_final_coding", "variant_classification"
    ))

    suppressWarnings(
      mutations_select <- mutations_select %>%
        mutate(consequence_final_coding = forcats::fct_recode(.data$variant_classification, !!!consequence_vec))
    )

    mutations_select_2 <- mutations_select %>%
      mutate(
        resource = resource, reference_genome = reference_genome,
        token = token,
        consequence_final_coding = tolower(.data$consequence_final_coding)
      ) %>%
      select(any_of(c(
        "resource",
        "hugo_symbol",
        "hgv_sp_short",
        "reference_genome",
        "consequence_final_coding",
        "protein_pos_start",
        "protein_pos_end",
        "token"
      ))) %>%
      mutate()


    requests <- mutations_select_2 %>%
      rowwise() %>%
      mutate(requests = list(oncokb_api_request(
        resource = resource,
        hugoSymbol = hugo_symbol,
        alteration = hgv_sp_short,
        referenceGenome = reference_genome,
        consequence = consequence_final_coding,
        proteinStart = protein_pos_start,
        proteinEnd = protein_pos_end,
        token = token
      ))) %>%
      pull(requests)

    list_of_responses <- map(requests, function(request) {
      req_perform(request) |>
        httr2::resp_body_string()
    })

    all_mutations_oncokb <- map_df(list_of_responses, parse_responses)
    all_mutations_oncokb <- bind_cols(
      mutations_select[, c("event_index", "sample_id")],
      all_mutations_oncokb
    )

    #**  Check errors ---------
    if ("message" %in% colnames(all_mutations_oncokb)) {
      messages <- stats::na.omit(all_mutations_oncokb$message)

      if (length(messages) > 0) {
        cli::cli_warn("{.val {length(messages)}} obervation(s) came back with the following information: {.code {unique(messages)}}
                        \nThese may be due to missing/NA data in columns used for annotation. To further check these values, see the {.col message} and {.col status}
                      columns of your resulting dataframe (e.g. {.code results[which(!is.na(results$message)), c('message', 'status')] })")
      }
    }
  }
  # Clean Results  ----------------------------------------------------------

  all_mutations_oncokb <- .clean_query_results(
    query_result = all_mutations_oncokb,
    return_simple = return_simple,
    return_query_params = return_query_params,
    original_data = mutations
  )

  # * Tumor Type - Remove Col if None  ------------------------------------------

  all_mutations_oncokb <- .tumor_type_warning(
    annotate_tumor_type = annotate_tumor_type,
    data = all_mutations_oncokb
  )

  return(all_mutations_oncokb)
}
