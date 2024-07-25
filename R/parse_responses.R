

parse_responses <- function(raw_response) {

    parsed_json <- jsonlite::fromJSON(raw_response,
                                      flatten = TRUE,
                                      simplifyVector = TRUE)

    parsed_tibble <- unlist(parsed_json, recursive=TRUE) %>%
      tibble::enframe() %>%
      tidyr::pivot_wider(names_from = "name",
                         values_fn = function(x) paste(x, collapse=","))

  }

