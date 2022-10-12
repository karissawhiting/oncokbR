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


