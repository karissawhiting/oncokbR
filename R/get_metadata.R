# Versions ----------------------------------------------------------------

get_oncokb_versions <- function() {

  # Set resource ----------------
  resource = "info"

  # Make request ---------------
  resp <- oncokbR_request(resource = resource)

  # Clean response --------------
  tibble::tibble(
    oncokb_api_version = resp$apiVersion$version,
    oncokb_data_version = resp$dataVersion$version,
    oncokb_app_version = resp$appVersion$version
  )
}



# Available Levels of Evidence --------------------------------------------

available_levels <- function() {

  # Set resource ----------------
  resource = "levels"

  # Make request ---------------
  resp <- oncokbR_request(resource = resource)

  # Clean response --------------
  tibble::tibble(
    oncokb_api_version = names(resp),
    oncokb_data_version = purrr::map_chr(resp, ~.x),
  )
}

