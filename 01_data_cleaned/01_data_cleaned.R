# Functions and packages ----
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

## names of packages we want installed (if not installed yet) and loaded
package_vec <- c(
  "rgbif",
  "tidyverse",
  "here"
)

## executing install & load for each package
sapply(package_vec, install.load.package)

# exploring datasets ----
load(here(
  "00_raw_data",
  "data_tetrapods_southamerica.RData")
)

# conferir as chaves dos dados baixados
# conferir % das colunas var. resp e preditoras
View(data_tetrapodstraits)
table(data_tetrapodstraits$occurrenceStatus)

View(data_tetrapodstraits %>% 
       filter(occurrenceStatus == " "))
View(data_tetrapodstraits %>% 
       filter(taxonRank == "UNRANKED")) 

# key of each occurence
anyDuplicated(data_tetrapodstraits$gbifID)

data_tetrapods_rgbif <- data_tetrapodstraits %>%
  select("speciesKey", "class", "order", "family", "species",
         "countryCode", "locality", "stateProvince",
         "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters",
         "coordinatePrecision", "eventDate", "year",
         "taxonRank", "occurrenceStatus", "gbifID") %>%
  #filter(taxonRank == "SPECIES" | taxonRank == "SUBSPECIES") %>%
  filter(taxonRank == "SPECIES") %>%
  filter(occurrenceStatus != "ABSENT") %>%
  select(-"occurrenceStatus", -"taxonRank") %>%
  mutate(across(everything(), replace_empty_with_na))


nrow(data_tetrapods_rgbif) # 1.965.555 occurences species
table(data_tetrapods_rgbif$class)
glimpse(data_tetrapods_rgbif)

# how much data is na?
data_tetrapods_rgbif %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year)), 
            missing_date = sum(is.na(eventDate)))

# nested dataframe
data_tetrapods_nested <- data_tetrapods_rgbif %>%
  group_by(species) %>%
  nest(event_table = c("gbifID", "eventDate", "year", "countryCode", "locality",
                       "stateProvince","decimalLatitude", "decimalLongitude",
                       "coordinateUncertaintyInMeters","coordinatePrecision")) %>%
  ungroup() %>%
  arrange("speciesKey", "class", "order", "family", "species")

anyDuplicated(data_tetrapods_nested$species) # nenhum nome repetido
View(data_tetrapods_nested %>% filter(species == "Boana polytaenia"))

View(data_tetrapods_nested %>%
     filter(duplicated(species)))
