
# Prep and Check Data --------------------------------------------------------
annotate_mutations <- function(mutations,
                               annotate_by = c("protein_change", "hgvsg"),
                               annotate_tumor_type = NULL,
                               return_simple = TRUE,
                               return_query_params = FALSE,
                               reference_genome = "GRCh37",
                               token = get_oncokb_token()
) {

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
         "protein_change" = c("sample_id",
                              "hugo_symbol",
                              "hgv_sp_short",
                              "variant_classification",
                              "protein_pos_start",
                              "protein_pos_end"))


  .check_required_cols(data = mutations,
                       required_cols = required_cols,
                       data_name = "mutations")


  # Check Parameters ----------------------------------------------------------

  # * `reference_genome`
  if (!inherits(reference_genome, "character")) {
    stop("`reference_genome` must be a character")
  }

  # Annotate Mutations --------------------------------------------------------

  # Set Resource
  resource <- switch(annotate_by,
                     "hgvsg" = "annotate/mutations/byHGVSg",
                     "protein_change" = "annotate/mutations/byProteinChange")

  # Check Token
  validate_oncokb_token(token = token)

  # create event index
  mutations <- mutate(mutations, event_index = 1:nrow(mutations))

  # subset data to only necessary columns
  mutations_select <- mutations %>%
    select(any_of(c("event_index", required_cols, "tumor_type")))

  # * If HGVSG ----------

  if(annotate_by == "hgvsg") {

    if (!inherits(mutations_select$hgv_sg, "character")) {
      stop("`hgv_sg` must be a character")}

    requests <- map(mutations_select$hgv_sg, ~oncokb_api_request(
      resource = resource,
      hgvsg = .x,
      referenceGenome = reference_genome,
      token = token
    ))

    list_of_responses <- map(requests, function(request) {
      req_perform(request) |>
        httr2::resp_body_string()})

    all_mutations_oncokb <- map_df(list_of_responses, parse_responses)
    all_mutations_oncokb <- bind_cols(mutations_select[, c('event_index', 'sample_id')],
                                      all_mutations_oncokb)

    #*  Check errors ---------
    messages <- na.omit(all_responses$message)

    if (length(messages) > 0) {
      cli::cli_warn("{.val {length(messages)}} obervation(s) came back with the following information: {.code {unique(messages)}}
                      \nThese may be due to missing/NA data in columns used for annotation. To further check these values, see the {.col message} and {.col status}
                    columns of your resulting dataframe (e.g. {.code results[which(!is.na(results$message)), c('message', 'status')] })")
    }

    return(all_responses)

   # LEFT OFFF HERE!!!!!
    # * If Protein Change -----

    if(annotate_by == "protein_change") {

      # Checks Specific to Protein Change -----------------

      # create a protein start/end column from protein position if needed
      mutations <- .check_protein_start_end(column_names, mutations)


      required_cols_pro <- c("sample_id",
                             "hugo_symbol",
                             "hgv_sp_short",
                             "variant_classification",
                             "protein_pos_start",
                             "protein_pos_end")

      .check_required_cols(data = mutations,
                           required_cols = required_cols_pro,
                           data_name = "mutations")

      # * Check Variant Consequence  -----------

      variant_options <- tolower(unique(stats::na.omit(unlist(oncokbR::consequence_map))))
      variant_in_data <- tolower(unique(mutations$variant_classification))

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

  # Clean Results  ----------------------------------------------------------

  all_mutations_oncokb <- .clean_query_results(
    query_result = all_mutations_oncokb,
    return_simple = return_simple,
    return_query_params = return_query_params,
    original_data = mutations)

  # * Tumor Type - Remove Col if None  ------------------------------------------

  all_mutations_oncokb <- .tumor_type_warning(
    annotate_tumor_type = annotate_tumor_type,
    data = all_mutations_oncokb)

  return(all_mutations_oncokb)
}


# Annotate Mutations ------------------------------------------------------


annotate_mutations_by_hgvsg <- function(
  hgvsg = c("7:g.140453136A>T"),
  referenceGenome = "GRCh37",
  annotate_tumor_type = NULL,
  token = get_oncokb_token()) {



  # Check parameters -----------------------------------------------------------


hgvsg = c("7:g.140453136A>T", NA, "3:g.179218294G>A")
referenceGenome = c("GRCh37", "GRCh37", "GRCh37")

x <- annotate_mutations_by_hgvsg(hgvsg = mutations$hgv_sg, referenceGenome = rep(c("GRCh37"), nrow(mutations)))
