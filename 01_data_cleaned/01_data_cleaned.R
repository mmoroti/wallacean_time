# Functions and packages ----
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

## names of packages we want installed (if not installed yet) and loaded
package_vec <- c(
  "cowplot",
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

# exists a difference between "year" and "event date" of 6041 occurences
# these differences is associate a occurences between imprecise dates
# check documentation https://encurtador.com.br/XmP8P
View(data_tetrapods_rgbif %>%
       filter(!is.na(eventDate) & is.na(year)))

# select occurences at
year_filter <- 2020

data_tetrapods_filter <- data_tetrapods_rgbif %>% 
  drop_na(decimalLatitude,
          decimalLongitude,
          year) %>%
  filter(year <= year_filter)

data_tetrapods_filter %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year)), 
            missing_date = sum(is.na(eventDate))) # without missing data =)

nrow(data_tetrapods_filter)

# comparative plots of original vs filtered dataframes 
df_filter <- as.data.frame(table(data_tetrapods_filter$year)) %>%
  rename(Year = Var1)
df_original <- as.data.frame(table(data_tetrapods_rgbif$year)) %>%
  rename(Year = Var1) 

# Convert the Year column to numeric if it's a factor
df_filter$Year <- as.numeric(as.character(df_filter$Year))
df_original$Year <- as.numeric(as.character(df_original$Year)) 

df_comparative <- full_join(df_original, df_filter,
          by = "Year") %>%
  rename(Original = Freq.x,  Filtered = Freq.y) %>%
  filter(Year <= year_filter) %>%
  pivot_longer(cols = c(Original, Filtered), 
               names_to = "Dataset",
               values_to = "Occurences") %>%
  replace_na(list(Occurences = 0))

ggplot(df_comparative, aes(x = Year, y = Occurences, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", width = 2, alpha = 0.6) +
  geom_line(stat = "smooth", aes(color = Dataset), size = 1.0, linetype = "dashed") + 
  theme_bw() +
  labs(title = "",
       y = "Frequency of occurences")

# nested dataframe
data_tetrapods_nested <- data_tetrapods_filter %>%
  group_by(species) %>%
  nest(event_table = c("gbifID", "eventDate", "year", "countryCode", "locality",
                       "stateProvince","decimalLatitude", "decimalLongitude",
                       "coordinateUncertaintyInMeters","coordinatePrecision")) %>%
  ungroup() %>%
  arrange(class, order, family, species)

# check data
table(data_tetrapods_nested$class)
table(data_tetrapods_nested$order)
anyDuplicated(data_tetrapods_nested$species) # nenhum nome repetido
anyDuplicated(data_tetrapods_nested$speciesKey)

# counting events and ordering years
data_tetrapods_nested$count_events <- NA

for (i in 1:nrow(data_tetrapods_nested)) {
  data_tetrapods_nested[i,"count_events"] <- nrow(
    data_tetrapods_nested$event_table[[i]])
  
  data_tetrapods_nested$event_table[[i]] <- 
    data_tetrapods_nested$event_table[[i]][
      order(data_tetrapods_nested$event_table[[i]]$year), 
    ]
}

data_tetrapods_nested %>%
  filter(count_events > 10)

save(data_tetrapods_nested, 
     file = here(
       "01_data_cleaned",
       "data_tetrapods_nested.RData")
)
