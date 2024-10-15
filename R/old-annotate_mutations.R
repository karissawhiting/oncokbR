#'
#'
#' #' Annotate Mutations
#' #'
#' #' Annotate maf and maf-like files using information about protein change or HGVSg.
#' #' You can annotated with tumor types indicated (oncotree codes), or without.
#' #'
#' #' @param mutations a mutations file in MAF, or similar format
#' #' @param annotate_by Can indicate whether to annotate by "protein_change" or "hgvsg (see oncoKB API docs for more info).
#' #' "Default is `protein_change`
#' #' @param return_simple Default is `TRUE` where only a set of the most common columns are returned from oncoKB annotator
#' #' see `oncokbR::output_dictionary` for more information on what is returned. If `FALSE` all raw columns are returned from API.
#' #' @param return_query_params If `TRUE`, the specific parameters used to query the API are returned in new columns.
#' #' This can be useful for troubleshooting the annotator. Default is `FALSE`.
#' #'
#' #' @return an annotated mutations file
#' #' @export
#' #' @import dplyr
#' #'
#' #' @examples
#' #' mutations <- oncokbR::blca_mutation[1:10, ]
#' #'
#' #' # Annotate without tumor type -----
#' #' x <- annotate_mutations(mutations)
#' #'
#' #' # Annotate with tumor type -----
#'
#' #'
#' #' # Annotate by HGVSg
#' #' z <- annotate_mutations(oncokbR::mutations_hgvsg, annotate_by =  "hgvsg")
#' #'
#' annotate_mutations <- function(mutations, annotate_by = c("protein_change", "hgvsg"),
#'                                return_simple = TRUE,
#'                                return_query_params = FALSE
#'                                ) {
#'
#'   # standardize column names
#'   mutations <- rename_columns(mutations)
#'   column_names <- colnames(mutations)
#'
#'   # check for tumor type
#'   annotate_tumor_type <- ("tumor_type" %in% column_names)
#'
#'   annotate_by <- match.arg(annotate_by)
#'
#'   # Protein Change -------------------------------------------------------------
#'
#'   if(annotate_by == "protein_change") {
#'
#'     # * Required Columns ------
#'
#'
#'     }
#'
#'
#'   # HGVSG ----------------------------------------------------------------------
#'
#'   # * Required Columns ------
#'
#'   if(annotate_by == "hgvsg") {
#'     required_cols_pro <- c("sample_id", "hgv_sg")
#'
#'     .check_required_cols(data = mutations,
#'                          required_cols = required_cols_pro,
#'                          data_name = "mutations")
#'
#'   }
#'
#'   # Annotate Mutations  -------------------------------------------------------
#'
#'
#'   mutations <- mutate(mutations, index = 1:nrow(mutations))
#'
#'   # all_mut_oncokb <- switch(annotate_by,
#'   #                          "protein_change" =
#'   #                            .annotate_mutations_by_protein_change(mutations = mutations,
#'   #                                                                  annotate_tumor_type = annotate_tumor_type),
#'   #                          "hgvsg" = .annotate_mutations_by_hgvsg(mutations = mutations,
#'   #                                                                 annotate_tumor_type = annotate_tumor_type)
#'   #                          )
#'
#'   if(annotate_by == "hgvsg") {
#'     all_mut_oncokb <- .annotate_mutations_by_hgvsg(hgvsg = mutations$hgv_sg,
#'                                 referenceGenome = rep("GRCh37", nrow(mutations)))
#'   }
#'   # Clean Results  ----------------------------------------------------------
#'
#'   all_mut_oncokb <- .clean_query_results(
#'     query_result = all_mut_oncokb,
#'     return_simple = return_simple,
#'     return_query_params = return_query_params,
#'     original_data = mutations)
#'
#'   # Tumor Type - Remove Cols if None  ------------------------------------------
#'
#'   all_mut_oncokb <- .tumor_type_warning(
#'     annotate_tumor_type = annotate_tumor_type,
#'     data = all_mut_oncokb)
#'
#'   return(all_mut_oncokb)
#'
#'
#' }
#'
#' # Sub Annotators By Data Type ----------------------------------------------------
#'
#' #' Annotate Mutations By Protein Change
#' #'
#' #' @param mutations
#' #'
#' #' @return a dataframe of annotated mutations
#' #' @export
#' #'
#' #' @keywords internal
#' #' @examples
#' #' mutations <- rename_columns(oncokbR::blca_mutation[1:10, ])
#' #' mutations$index <- 1:nrow(mutations)
#' #' mutations$consequence_final_coding = "any"
#' #' .annotate_mutations_by_protein_change(mutations, annotate_tumor_type = FALSE)
#' #'
#' .annotate_mutations_by_protein_change <- function(mutations, annotate_tumor_type) {
#'
#'   all_mut_oncokb <- mutations %>%
#'     select(any_of(c("index", "sample_id", "hugo_symbol", "hgv_sp_short", "consequence_final_coding",
#'                     "protein_pos_start", "protein_pos_end", "tumor_type")))
#'
#'   make_url <- function(index, sample_id,
#'                        hugo_symbol,
#'                        hgv_sp_short,
#'                        consequence_final_coding,
#'                        protein_pos_start, protein_pos_end, tumor_type) {
#'
#'     url <- glue::glue("https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?hugoSymbol=",
#'                       "{hugo_symbol}&alteration={hgv_sp_short}&referenceGenome=GRCh37&consequence=",
#'                       "{consequence_final_coding}&proteinStart={protein_pos_start}&proteinEnd={protein_pos_end}")
#'
#'
#'     if(annotate_tumor_type) {
#'       url <- glue::glue(url, "&tumorType={tumor_type}")
#'     }
#'
#'     token <- Sys.getenv('ONCOKB_TOKEN')
#'     resp <- httr::GET(url,
#'                       httr::add_headers(Authorization = paste("Bearer ",
#'                                                               token,
#'                                                               sep = "")))
#'
#'     parsed <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
#'                                  flatten = TRUE,
#'                                  simplifyVector = TRUE)
#'
#'     # parsed <- parsed %>% discard("query")
#'     parsed <-unlist(parsed, recursive=TRUE) %>%
#'       tibble::enframe() %>%
#'       tidyr::pivot_wider(names_from = "name",
#'                          values_fn = function(x) paste(x, collapse=","))
#'
#'     parsed$sample_id <- sample_id
#'     parsed$index <- index
#'
#'     parsed
#'
#'   }
#'
#'   purrr::pmap_df(all_mut_oncokb,  make_url)
#'
#'
#'  }
#'
#'
#'
#' #' Annotate Mutations By HGVSg
#' #'
#' #' @param mutations
#' #'
#' #' @return a dataframe of annotated mutations
#' #' @export
#' #'
#' #' @keywords internal
#' #' @examples
#' #'
#' #'
#' #' mutations <- rename_columns(oncokbR::mutations_hgvsg)
#' #' mutations$index = 1:nrow(mutations)
#' #' .annotate_mutations_by_hgvsg(mutations, annotate_tumor_type = FALSE)
#' #'
#' .annotate_mutations_by_hgvsg <- function(mutations, annotate_tumor_type) {
#'
#'   all_mut_oncokb <- mutations %>%
#'     select(any_of(c("index", "sample_id", "hgv_sg", "tumor_type")))
#'
#'   make_url <- function(index, sample_id, hgv_sg, tumor_type) {
#'
#'     url <- glue::glue("https://www.oncokb.org/api/v1/annotate/mutations/byHGVSg?hgvsg=",
#'                       "{hgv_sg}&referenceGenome=GRCh37")
#'
#'     if(annotate_tumor_type) {
#'       url <- glue::glue(url, "&tumorType={tumor_type}")
#'     }
#'
#'     url <- utils::URLencode(url)
#'
#'     token <- Sys.getenv('ONCOKB_TOKEN')
#'     resp <- httr::GET(url,
#'                       httr::add_headers(Authorization = paste("Bearer ",
#'                                                               token,
#'                                                               sep = "")))
#'
#'     parsed <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
#'                                  flatten = TRUE,
#'                                  simplifyVector = TRUE)
#'     # parsed <- parsed %>% discard("query")
#'     parsed <- unlist(parsed, recursive=TRUE) %>%
#'       tibble::enframe() %>%
#'       tidyr::pivot_wider(names_from = "name",
#'                          values_fn = function(x) paste(x, collapse=","))
#'
#'     parsed$sample_id <- sample_id
#'     parsed$index <- index
#'     parsed
#'
#'   }
#'
#'   purrr::pmap_df(all_mut_oncokb,  make_url)
#'
#'
#' }
#'
#'
#'
