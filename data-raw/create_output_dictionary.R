library(dplyr)
library(readr)

output_dictionary <- read_csv(
  here::here("data-raw", "output-data-dictionary.csv"),
  trim_ws = TRUE)


usethis::use_data(output_dictionary, overwrite = TRUE)
