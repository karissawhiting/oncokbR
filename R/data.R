
#' Mutation data from Nonmuscle invasive bladder cancer patients (MSK Eur Urol 2017)
#'
#' A dataset of mutations, sourced from cBioPortal
#'
#' @format A data frame with 1595 rows of mutation data from

"blca_mutation"

#' Copy number alteration data from Nonmuscle invasive bladder cancer patients (MSK Eur Urol 2017)
#'
#' A dataset of CNA, sourced from cBioPortal
#'
#' @format A data frame with 168 rows of CNA data from
"blca_cna"


#' Mutation data from Nonmuscle invasive bladder cancer patients (MSK Eur Urol 2017)
#'
#' A dataset of mutations, sourced from cBioPortal
#'
#' @format A data frame with 13 rows of structural variant data from
"blca_sv"

#' Generated mutation data
#'
#' A dataset of mutations, with hgvsg data
#'
#' @format A data frame with 13 rows of structural variant data from
"mutations_hgvsg"

#' Mutation data from ACC TCGA
#'
#' A dataset of mutations, sourced from cBioPortal
#'
#' @format A data frame with 13486 rows of mutation data
"acc_mutation"

#' Copy number alteration data from ACC TCGA
#'
#' A dataset of CNA, sourced from cBioPortal
#'
#' @format A data frame with 23878 rows of CNA data from
"acc_cna"

#' Consequence Map
#'
#' Dataset of recoding values for consequence mutation data
#'
#' @format A data frame with 23878 rows of CNA data from
"consequence_map"

#' Data Frame of Column Names
#'
#' Data frame of accepted data names for standard genomic files. This serves as a
#' dictionary to help disambiguate raw column names from user entered mutation,
#' CNA or structural variant data
#'
#' @format A data frame
#' \describe{
#'     \item{maf_column_name}{data field names as they appear in common MAF file}
#'     \item{api_column_name}{data field names as they appear in common cBioPortal API retrieved files}
#'     \item{mutation_input}{does this field appear in mutation files?}
#'     \item{fusion_input}{does this field appear in mutation/sv files?}
#'     \item{cna_input}{does this field appear in CNA files?}
#'     \item{definition}{variable definition}
#'     \item{notes}{data notes}
#'     \item{sc_maf_column_name}{snake case version of `maf_column_name`}
#'     \item{sc_api_column_name}{snake case version of `api_column_name`}
#'     \item{caps_maf_column_name}{upper case version of `maf_column_name`}
#'     \item{caps_api_column_name}{upper case version of `api_column_name`}
#'     \item{internal_column_name}{name used for each field for all internal processing functions}
#'     }
"names_df"
