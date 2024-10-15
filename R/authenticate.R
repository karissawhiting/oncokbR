

# Authentication -----------------------------------------------------------

#' Get oncoKB Access Token
#'
#' Convenience function that retrieves oncoKB token System Environment variable "ONCOKB_TOKEN"
#' @export
#' @return Returns a string with token if successfully authenticated, or
#' a warning that token could not be found.
#' @author Karissa Whiting
#' @examples
#' \dontrun{
#' get_oncokb_token()
#' }
#'
get_oncokb_token <- function(token = NULL) {

  token <- token %||% Sys.getenv("ONCOKB_TOKEN")

  if (identical(token, "")) {
    rlang::warn("No ONCOKB_TOKEN in `.Renviron`. Try `usethis::edit_r_environ()` to add a token")
  }
  token
}


# Check Token -------------------------------------------------------------

validate_oncokb_token <- function(token = get_oncokb_token()) {

  resource <- c("tokens", token)

  query <- oncokb_api_request(
    resource = resource,
    api_url = "https://www.oncokb.org/api/",
    token = get_oncokb_token()
  )

  # Custom error ---may want to make this arg in main API function
  error_body <- function(resp) {
    paste0("There was an error validating your token. Details: ",
    resp_body_json(resp)$error,
    resp_body_json(resp)$detail,
    " Go to `https://www.oncokb.org/` to get a valid token.")
  }

  response <- query |>
    req_error(body = error_body) |>
    req_perform(verbosity = 0) |>
    resp_body_json()

  expiration_date <- as.Date(response$expiration)

  if(expiration_date < as.Date(Sys.Date())) {
    cli::cli_alert_danger("Token has expired. Go to https://www.oncokb.org/ to renew. ")
  }

  cli::cli_alert_success(
    "Token valid for {.val {response$user$login}}")
  cli::cli_alert_info(
    "Token expiration date: {.val {expiration_date}}")
}

