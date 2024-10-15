#' Rename columns from API results to work with gnomeR functions
#'
#' @param df_to_check a data frame for which to check colums
#'
#' @return a renamed data frame
#' @export
#' @examples
#'
#' rename_columns(df_to_check = oncokbR::blca_mutation)
#' rename_columns(df_to_check = oncokbR::blca_sv)
#'
rename_columns <- function(df_to_check) {

  names_df_long <- oncokbR::names_df %>%
    select(contains("_column_name")) %>%
    tidyr::pivot_longer(-"internal_column_name") %>%
    distinct()

  which_to_replace <- intersect(names(df_to_check), unique(names_df_long$value))

  # create a temporary dictionary as a named vector- this should have all relevant values, including those unchanged
  names_dict <- names_df_long %>%
    dplyr::filter(.data$value %in% which_to_replace) %>%
    select("internal_column_name",  "value") %>%
    dplyr::distinct() %>%
    tibble::deframe()


  if(length(names_dict) > 0) {

    # store details on what has been changed.
    message <- purrr::map2_chr(names(names_dict),
                               names_dict,
                               ~paste0(.y, " renamed ", .x))

    names(message) <- rep("!", times = length(message))


    # rename those variables only
    df_to_check <- df_to_check %>%
      dplyr::rename(!!names_dict)

    attr(df_to_check, "names_dict") <- names_dict
  }

  return(df_to_check)
}


#' Internal function to recode numeric CNA alteration values to factor values
#'
#' @param alteration_vector a vector of CNA alterations coded with any of the
#' following levels: neutral, deletion, amplification, gain, loss, homozygous deletion,
#' hemizygous deletion, loh, gain, high level amplification, 0, -1, -1.5, -2, 1, 2.
#'
#' @return a recoded CNA data set with factor alteration values. See details for code dictionary
#'
#' @details
#'
#' CNA is coded to the following key based on key: values below
#' - "neutral":  "0", "neutral",
#' - "deletion": "homozygous deletion", "-2",
#' - "deletion": "loh", "-1.5",
#' - "deletion": "hemizygous deletion", "-1",
#' - "amplification": "gain", "1",
#' - "amplification": high level amplification", "2",
#' @export
#' @examples
#' recode_cna(blca_cna$alteration[1:10])

recode_cna <- function(alteration_vector){

  # *TODO - should we auto-recode Unknown/Unk to NA? - need to think on it

  # General Checks -------------------------------------------------------------
  # *TODO (I may remove these if this is just an internal function TBD )
  # (as it's already done in higher level function) but may change

  alteration_vector = as.character(alteration_vector)

  # CNA Levels Checks -----------------------------------------------------------
  levels_in_data <- tolower(names(table(alteration_vector)))

  # source: https://docs.cbioportal.org/file-formats/#data-file-1
  # python annotator ref with codes: https://github.com/oncokb/oncokb-annotator/blob/47e4a158ee843ead75445982532eb149db7f3106/AnnotatorCore.py#L158
  allowed_cna_levels <- tibble::tribble(
    ~detailed_coding, ~numeric_coding,   ~final_coding,
    "neutral",          "0",        "neutral",
    "deep loss",          "-2",      "deletion",
    "deep loss",          "-1.5",    "deletion",
    "hemizygous deletion",          "-1",       "loss",
    "gain",          "1",        "gain",
    "high level amplification",          "2",   "amplification")



  all_allowed <- unlist(allowed_cna_levels)
  not_allowed <- levels_in_data[!levels_in_data %in% all_allowed]

  if(length(not_allowed) > 0) {
    cli::cli_abort(c("Unknown values in {.field alteration} field: {.val {not_allowed}}",
                     "Must be one of the following: {.val {all_allowed}}"))
  }

  # Recode CNA Levels  ----------------------------------------------------------

  # create a named vector for recoding to final coding
  recode_values <- c(allowed_cna_levels$detailed_coding, allowed_cna_levels$numeric_coding)
  names(recode_values) <- c(allowed_cna_levels$final_coding, allowed_cna_levels$final_coding)

  recoded_alterations <- suppressWarnings(
    forcats::fct_recode(alteration_vector, !!!recode_values)
  )


  return(recoded_alterations)
}



#' Helper to select columns for output of annotation functions
#'
#' @param results input annotation result data frame
#' @inheritParams annotate_mutations
#'
#' @return a string of select column names
#' @export
#'
#' @examples
.get_select_columns <- function(results,
                        return_simple,
                        return_query_params) {

  all_colnames <- names(results)
  query_cols <- all_colnames[stringr::str_detect(all_colnames, "oncokb_query_")]

  simple_cols <- switch(return_simple + 1,
                        all_colnames[!(all_colnames %in% query_cols)],
                        oncokbR::output_dictionary$output_column_name[which(oncokbR::output_dictionary$include_in_simple_output == "yes")])

  query_cols <- switch(return_query_params + 1,
                       NULL,
                       query_cols)

  selected_columns <- unique(c(query_cols, simple_cols))

  return(selected_columns)

}



#' Clean API Query Results
#'
#' @param query_result a data frame returned from API
#' @param original_data the original user input data frame
#' @inheritParams annotate_mutations
#' @return a cleaned data frame
#' @export
#'
.clean_query_results <- function(query_result,
                                 original_data,
                                 return_simple,
                                 return_query_params) {

  query_result <- query_result %>%
    janitor::clean_names() %>%
    dplyr:: rename_with(~paste0("oncokb_", .),
                        .cols = -c("sample_id", "event_index"))

  select_oputput_columns <- .get_select_columns(query_result,
                                               return_simple = return_simple,
                                               return_query_params = return_query_params)
  query_result <- select(query_result,
                           any_of(c("sample_id", "event_index",
                                    select_oputput_columns)))

  query_result <- query_result %>%
    left_join(original_data, ., by = c("sample_id", "event_index"))

  return(query_result)
}


#' Check a Data Frame for Required Columns
#'
#' @param data A data frame to check
#' @param required_cols A character specifying names of columns to check
#' @param data_name Optionally specify how the data set should be called in error message.
#' Default is NULL and will call it a generic name.
#' @return If data set doesn't have required columns it will return an error message.
#' If it does have required columns, nothing will be returned
#' @keywords internal

.check_required_cols <- function(data, required_cols, data_name = NULL) {

  data_name <- data_name %||% ""
  column_names <- colnames(data)
  which_missing <- required_cols[which(!(required_cols %in% column_names))]

  if(length(which_missing) > 0) {
    cli::cli_abort("The following required columns are missing in your {data_name} data: {.field {which_missing}}")
  }

}


#' Warn if no tumor type annotation possible
#'
#' @param annotate_tumor_type TRUE if tumor type column is available in data
#' @param data mutation, cna or sv data
#'
#' @return a warning and data set with edited columns
#' @export
#'
#' @examples
.tumor_type_warning <- function(annotate_tumor_type,
                                data) {

  if (!annotate_tumor_type) {

    data <- data %>%
      select(-contains("treatments"))

    cli::cli_alert_info("No {.val tumor_type} found in data. No treatment-level annotations will be returned.")
  }

  return(data)

}





# Create Protein Position Column -------------------------------------------
.check_protein_start_end <- function(column_names, mutations) {

  if(!("protein_pos_start" %in% column_names | "protein_pos_end" %in% column_names) &
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

    return(mutations)
}



# Takes a vector as input
.check_consequence <- function(variant_classification) {

  # * Check Variant Consequence  -----------

  variant_options <- tolower(unique(stats::na.omit(unlist(oncokbR::consequence_map))))
  variant_in_data <- tolower(unique(variant_classification))

#  not_allowed <- stats::na.omit(variant_in_data[!(variant_in_data %in% variant_options)])
  not_allowed <- stats::na.omit(setdiff(variant_in_data, variant_options))

  # Maybe turn into warning
  if(length(not_allowed) > 0) {
    cli::cli_abort("The following variant classification levels are not recognized: {.code {not_allowed}}.
                     Please remove or recode these to continue (see {.code oncokbR::consequence_map} for allowed values)")
  }
}
