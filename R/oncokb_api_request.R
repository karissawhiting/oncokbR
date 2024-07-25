


oncokb_api_request <- function(resource = "info",
                               api_url = "https://www.oncokb.org/api/v1",
                               ...,
                               token = NULL) {

  params <- list(...)

  request(api_url) |>
    req_url_path_append(resource) |>
    req_url_query(!!!params) |>
    req_auth_bearer_token(token) |>
    req_error(is_error = \(resp) FALSE)
}


