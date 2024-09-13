# TODO adaptar para usar o ambiente virtual
# Functions and packages ----
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    renv::restore()
    #install.packages(x, repos = "http://cran.us.r-project.org")
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
rm(list=ls()); gc() # clean local enviroment

# load data
load(file.path(
  "00_raw_data",
  "tetrapodstraits_data.RData")
) 
load(file.path(
  "00_raw_data",
  "biotime_data.RData")
)
load(file.path(
  "00_raw_data",
  "splinks_data.RData")
)

# GBIF cleaned ----
# conferir as chaves dos dados baixados
# conferir % das colunas var. resp e preditoras
View(data_tetrapodstraits)

table(data_tetrapodstraits$occurrenceStatus) # remover ausencias

table(data_tetrapodstraits$taxonRank)
View(data_tetrapodstraits %>% 
       filter(taxonRank == "UNRANKED"))

# key of each occurence
anyDuplicated(data_tetrapodstraits$gbifID) 

data_tetrapods_clean_v1 <- data_tetrapodstraits %>%
  select("speciesKey", "class", "order", "family", "species",
         "countryCode", "locality", "stateProvince",
         "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters",
         "coordinatePrecision", "eventDate", "year",
         "taxonRank", "occurrenceStatus", "gbifID") %>%
  filter(taxonRank == "SPECIES" | taxonRank == "SUBSPECIES") %>%
  filter(occurrenceStatus != "ABSENT") %>%
  select(-"occurrenceStatus", -"taxonRank") %>%
  mutate(origin_of_data = "gbif")

nrow(data_tetrapods_clean_v1) # 2.489.856 occurences species
table(data_tetrapods_clean_v1$class)

#View(data_tetrapods_clean_v1 %>% 
#       filter(class == "")) # todas que nao tem classe correspondem ao mesmo sp

data_tetrapods_clean <- data_tetrapods_clean_v1 %>%
  mutate(class = if_else(class == "", "Squamata", class)) %>%
  mutate(eventDate = na_if(eventDate, ""))

table(data_tetrapods_clean$class)
glimpse(data_tetrapods_clean)
names(data_tetrapods_clean)

# how much data is na?
data_tetrapods_clean %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year)), 
            missing_date = sum(is.na(eventDate)))

# exists a difference between "year" and "event date" of 7280 occurences
# these differences is associate a occurences between imprecise dates
# check documentation https://encurtador.com.br/XmP8P
View(data_tetrapods_clean %>%
       filter(!is.na(eventDate) & is.na(year)))

data_tetrapods_filter <- data_tetrapods_clean %>% 
  drop_na(decimalLatitude,
          decimalLongitude,
          year) 

data_tetrapods_filter %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year)), 
            missing_date = sum(is.na(eventDate))) # without missing data =)

nrow(data_tetrapods_filter)
length(unique(data_tetrapods_filter$species))

# comparative plots of original vs filtered dataframes 
df_filter <- as.data.frame(table(data_tetrapods_filter$year)) %>%
  rename(Year = Var1)
df_original <- as.data.frame(table(data_tetrapods_clean$year)) %>%
  rename(Year = Var1) 

# Convert the Year column to numeric if it's a factor
df_filter$Year <- as.numeric(as.character(df_filter$Year))
df_original$Year <- as.numeric(as.character(df_original$Year)) 

df_comparative <- full_join(df_original, df_filter,
          by = "Year") %>%
  rename(Original = Freq.x,  Filtered = Freq.y) %>%
  pivot_longer(cols = c(Original, Filtered), 
               names_to = "Dataset",
               values_to = "Occurences") %>%
  replace_na(list(Occurences = 0))

ggplot(df_comparative, aes(x = Year, y = Occurences, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", width = 2, alpha = 0.6) +
  geom_line(stat = "smooth", aes(color = Dataset), size = 1.0, linetype = "dashed") + 
  theme_bw() +
  labs(title = "GBIF Dataset",
       y = "Frequency of occurences")

# Biotime cleaned ----
data_biotime_clean_v1 <- biotime_data_key_precleaned %>%
  select("speciesKey", "class", "order", "family", "species", 
         "LATITUDE", "LONGITUDE", "YEAR", "STUDY_ID") %>%
  rename(decimalLatitude = LATITUDE, 
         decimalLongitude = LONGITUDE,
         year = YEAR) %>%
  mutate(origin_of_data = "biotime")

names(data_biotime_clean_v1)

data_biotime_clean_v1 %>% 
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year))) # without missing data

# splink cleaned ----
names(splink_data_key_precleaned)

data_splink_clean_v1 <- splink_data_key_precleaned %>%
  select("speciesKey", "taxonclass", "ordem", "family", "scientificname", 
         "latitude", "longitude", "yearcollected") %>%
  rename(decimalLatitude = latitude, 
         decimalLongitude = longitude,
         class = taxonclass,
         order = ordem,
         year = yearcollected) %>%
  mutate(origin_of_data = "splink") %>%
  mutate(year = na_if(year, ""))
  
names(data_splink_clean_v1)
nrow(data_splink_clean_v1) # 660.781 ocorrencias

# alguns casos tem 0 em uma das coordenadas
# latitude com 0
nrow(data_splink_clean_v1 %>%
       filter(decimalLatitude == 0 & decimalLongitude != 0)) 
# longitude com 0
nrow(data_splink_clean_v1 %>%
       filter(decimalLatitude != 0 & decimalLongitude == 0)) 

nrow(data_splink_clean_v1 %>%
  filter(year == 0)) # nenhum ano zero

nrow(data_splink_clean_v1 %>%
  filter(decimalLatitude == 0 & decimalLongitude == 0)) # coordenadas com zero
# ocorrencias de coordenadas e longitudes iguais a zero = 243.006

nrow(data_splink_clean_v1 %>% 
  filter(is.na(year))) # 63.167 anos com NA

# remove coordenadas com 0 na lat e na long
# filtra dados sem ano
data_splink_clean <- data_splink_clean_v1 %>%
  filter(!(decimalLatitude == 0 & decimalLongitude == 0)) %>% 
  filter(!is.na(year))

nrow(data_splink_clean_v1) - nrow(data_splink_clean) # 268742 occ removed

data_splink_clean %>% 
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year))) # without missing data

# unir os dataframes ---- 
# a maioria sao de especies fora da america do sul
# por isso estamos perdendo tantos dados aqui
# quando for para fazer para o mundo sera necessario reconferir
#key_losing <- left_join(data.frame(speciesKey = chaves), 
#                        data_biotime_clean_v1,
#                        by = "speciesKey")
#
#length(unique(data_tetrapods_filter$speciesKey))
#length(unique(data_biotime_clean_v1$speciesKey))
#length(chaves)
#nrow(key_losing) # 1022095 occurences
#View(key_losing %>%
#  distinct(speciesKey, .keep_all = TRUE))
chaves_perdidas_biotime <- setdiff(
  data_biotime_clean_v1$speciesKey,
  data_tetrapods_filter$speciesKey
) 

data_biotime_clean_sa <- data_biotime_clean_v1 %>% 
  filter(!(speciesKey %in% chaves_perdidas_biotime))

chaves_perdidas_splink <- setdiff(
  data_splink_clean$speciesKey,
  data_tetrapods_filter$speciesKey
) 

data_splink_clean_sa <- data_splink_clean %>% 
  filter(!(speciesKey %in% chaves_perdidas_splink))

data_wallacean_knownledge <- bind_rows(
  data_tetrapods_filter,
  data_biotime_clean_sa,
  data_splink_clean_sa)

# TODO
# executar o coordsclean e limpar coordenadas erradas
# nested_df? 
# salvar arquivos para proximo passo 02_data_filter_polygons

# TODO needs adjusments ---
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

# Tue Sep 10 13:25:55 2024 ------------------------------
# TODO add biotime dataset
biotime_data_nested <- biotime_data_key %>%
  group_by(speciesKey) %>%
  nest(biotime_table = c("STUDY_ID", "LATITUDE", "LONGITUDE", "YEAR")) %>%
  ungroup() 

View(biotime_data_nested)

teste <- full_join(data_tetrapods_nested,
                   biotime_data_nested[,-1], 
                   by="speciesKey")

View(teste)


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
