

annotate_mutations <- function(mutations) {

  mutations <- rename_columns(mutations)

  all_mut_oncokb <- mutations %>%
    select(hugo_symbol, proteinChange, consequence_final_coding,
           proteinPosStart, proteinPosEnd) %>%
    mutate(tumorType = "BLCA")

  make_url <- function(hugoGeneSymbol, proteinChange,
                       consequence_final_coding,
                       proteinPosStart, proteinPosEnd, tumorType) {
    url <- glue::glue("https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?hugoSymbol=",
                      "{hugoGeneSymbol}&alteration={proteinChange}&referenceGenome=GRCh37&consequence=",
                      "{consequence_final_coding}&proteinStart={proteinPosStart}&proteinEnd={proteinPosEnd}&tumorType={tumorType}")

    token <- Sys.getenv('ONCOKB_TOKEN')
    resp <- httr::GET(url,
                      httr::add_headers(Authorization = paste("Bearer ",
                                                              token,
                                                              sep = "")))

    parsed <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                                 flatten = TRUE,
                                 simplifyVector = TRUE)
    # parsed <- parsed %>% discard("query")
    parsed <-unlist(parsed, recursive=TRUE) %>%
      enframe() %>%
      pivot_wider(names_from = name,
                  values_fn = function(x) paste(x, collapse=","))

    parsed

  }

  # Error in `dplyr::bind_rows()`:
  #   ! Can't combine `..1$query.referenceGenome` <character> and `..27$query.referenceGenome` <list>.
  #
  c <- pmap_df(all_mut_oncokb,  make_url)


  c <- c %>%
    janitor::clean_names()

  c2 <- c%>%
    select(-contains("query."))

  c3 <-c2 %>% select(query_hugo_symbol, query_alteration, query_consequence, gene_exist, variant_exist, oncogenic, mutation_effect_known_effect,
                     mutation_effect_description, hotspot,
                     gene_summary, variant_summary,
                     highest_sensitive_level, treatments_level)


  c3$oncogenic %>% table()

  c4 <- c3
  #%>% filter(oncogenic %in% c("Likely Oncogenic", "Oncogenic"))


  table(c4$highest_sensitive_level)

  c4 %>% select(highest_sensitive_level, treatments_level) %>%
    tbl_summary()

  all_mut2 <- all_mut %>%
    bind_cols(c4)

  all_mut2 <- all_mut2 %>%
    filter(oncogenic %in% c("Likely Oncogenic", "Oncogenic"))

  all_mut2 %>% select(highest_sensitive_level, treatments_level) %>%
    tbl_summary()
}
