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
    tidyr::pivot_longer(-.data$internal_column_name)


  which_to_replace <- intersect(names(df_to_check), unique(names_df_long$value))

  # create a temporary dictionary as a named vector
  temp_dict <- names_df_long %>%
    dplyr::filter(.data$value %in% which_to_replace) %>%
    select("internal_column_name",  "value") %>%
    dplyr::distinct() %>%
    tibble::deframe()


  if(length(temp_dict) > 0) {

    # store details on what has been changed.
    message <- purrr::map2_chr(names(temp_dict),
                               temp_dict,
                               ~paste0(.y, " renamed ", .x))

    names(message) <- rep("!", times = length(message))


    # rename those variables only
    df_to_check %>%
      dplyr::rename(!!temp_dict)
  }
}

#' Internal function to recode numeric CNA alteration values to factor values
#'
#' @param cna a maf (long) form data set of CNAs. Must include an alteration column.
#'
#' @return a recoded CNA data set with factor alteration values
#'
#'


.recode_cna_alterations <- function(cna){


  if(!("alteration" %in% colnames(cna))) {
    cli::cli_abort("An alteration column is missing from your cna data. Use pivot_cna_longer(),
                   which will recode alterations, instead if dataset is in API format.")
  }

  # Make sure hugo & alteration is character
  cna <- cna %>%
    mutate(hugo_symbol = as.character(.data$hugo_symbol)) %>%
    mutate(alteration = tolower(str_trim(as.character(.data$alteration))))


  #assess levels of alteration
  levels_in_data <- names(table(cna$alteration))

  allowed_chr_levels <- c(
    "neutral" = "0",
    "homozygous deletion" = "-2",
    "loh" = "-1.5", #this is a placeholder until methods cleared up
    "hemizygous deletion" = "-1",
    "gain" = "1",
    "high level amplification" = "2"
  )

  #pull any numbers not in allowed list out
  all_allowed <- c(allowed_chr_levels, names(allowed_chr_levels))
  not_allowed <- levels_in_data[!levels_in_data %in% all_allowed]

  #abort if unknown values exist
  if(length(not_allowed) > 0) {
    cli::cli_abort(c("Unknown values in {.field alteration} field: {.val {not_allowed}}",
                     "Must be one of the following: {.val {all_allowed}}"))
  }

  # recode the alteration varaible as factor with those levels
  # and suppress warnings on this
  suppressWarnings(
    cna <- cna %>%
      mutate(alteration = forcats::fct_recode(.data$alteration, !!!allowed_chr_levels))
  )

  return(cna)
}
