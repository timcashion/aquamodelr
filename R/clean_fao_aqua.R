#' Assembles and Cleans FAO's Global Fishery and Aquaculture Production Data from scratch
#'
#' \code{clean_fao_aqua}
#' Add on to rebuild_fish which 'Extracts data from FAO zipped file and merges them into "tidy" (i.e., long) format'
#' @param path_to_zipfile Name of zipped file of FAO data (if in current directory); otherwise, specify path to file
#' @return A merged, tidy dataset.
#' @examples
#' clean_fao_aqua("./raw_data/GlobalProduction_2019.1.0.zip")

#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom readr write_csv
#' @importFrom magrittr `%>%`
#' @export

clean_fao_aqua <- function(path_to_zipfile=NA, output_path=NA){
  source("./functions/rebuild_fish.R")
  fao_raw <- rebuild_fish(path_to_zipfile)
  cols <- c("country_name_en", "species", "year", "quantity", "source_code", "species_name_en", "species_scientific_name", "isscaap_group")
  aqua <- fao_raw %>%
    select(cols) %>%
    filter(source_code != "CAPTURE") %>%
    filter(year > 1994)
  readr::write_csv(aqua, )
}
