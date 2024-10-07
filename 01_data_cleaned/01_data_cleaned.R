# Functions and packages ----
comparative_plot <- function(data_original, data_filtered,
                             year_column, title = "Comparative Plot") {
  
  # Criar tabelas de frequência para os anos
  df_filter <- as.data.frame(table(data_filtered[[year_column]])) %>%
    rename(Year = Var1)
  df_original <- as.data.frame(table(data_original[[year_column]])) %>%
    rename(Year = Var1)
  
  # Converter a coluna 'Year' para numérico, caso seja fator
  df_filter$Year <- as.numeric(as.character(df_filter$Year))
  df_original$Year <- as.numeric(as.character(df_original$Year))
  
  # Fazer o join e preparar o dataframe comparativo
  df_comparative <- full_join(df_original, df_filter, by = "Year") %>%
    rename(Original = Freq.x, Filtered = Freq.y) %>%
    pivot_longer(cols = c(Original, Filtered), names_to = "Dataset", values_to = "Occurences") %>%
    replace_na(list(Occurences = 0))
  
  # Criar o gráfico comparativo
  ggplot(df_comparative, aes(x = Year, y = Occurences, fill = Dataset)) +
    geom_bar(stat = "identity", position = "dodge", width = 2, alpha = 0.6) +
    geom_line(stat = "smooth", aes(color = Dataset), linewidth = 1.0, linetype = "dashed") + 
    theme_bw() +
    labs(title = title, y = "Frequency of occurrences")
}

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
  "here",
  "CoordinateCleaner",
  "sf",
  "countrycode"
)

# executing install & load for each package
# TODO adaptar para uso do ambiente virtual renv
sapply(package_vec, install.load.package)

# GBIF FILTER AND CLEAN ----
rm(list=ls()); gc() # clean local enviroment
# load data
start_time <- Sys.time()
load(file.path(
  "00_raw_data",
  "tetrapodstraits_data.RData")
)

# conferir as chaves dos dados baixados
# conferir % das colunas var. resp e preditoras
#head(data_tetrapodstraits)
#table(data_tetrapodstraits$occurrenceStatus) # remover ausencias
#table(data_tetrapodstraits$taxonRank)

# key of each occurence
anyDuplicated(data_tetrapodstraits$gbifID) 
nrow(data_tetrapodstraits) # 2513477 occurences species
table(data_tetrapodstraits$class)

data_tetrapods_clean <- data_tetrapodstraits %>%
  filter(taxonRank == "SPECIES" | taxonRank == "SUBSPECIES") %>%
  filter(occurrenceStatus != "ABSENT") %>%
  mutate(origin_of_data = "gbif") %>%
  mutate(eventDate = na_if(eventDate, "")) %>%
  mutate(class = if_else(class == "", "Squamata", class)) # TODO REVER quando usar dados globais

nrow(data_tetrapods_clean) # 2.489.856 occurences species
table(data_tetrapods_clean$class)

# how much data is na?
data_tetrapods_clean %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year)), 
            missing_date = sum(is.na(eventDate)))

# exists a difference between "year" and "event date" of 7280 occurences
# these differences is associate a occurences between imprecise dates
# check documentation https://encurtador.com.br/XmP8P
#View(data_tetrapods_clean %>%
#       filter(!is.na(eventDate) & is.na(year)))
data_tetrapods_filter <- data_tetrapods_clean %>% 
  drop_na(decimalLatitude,
          decimalLongitude,
          year) %>%
  select(-eventDate)

data_tetrapods_filter %>%
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year))) # without missing data =)
            #missing_date = sum(is.na(eventDate))) 

nrow(data_tetrapods_filter) # 1153622 total occurences 
table(data_tetrapods_filter$class) # occurences per classes 

# coordsclean
# convert country code from ISO2c to ISO3c
data_tetrapods_filter$countryCode <-  countrycode(
  data_tetrapods_filter$countryCode,
  origin =  'iso2c',
  destination = 'iso3c')

flags <- clean_coordinates(x = data_tetrapods_filter,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals", "centroids", "equal",
                            "gbif","institutions", "outliers",
                            "seas","zeros")) # Flagged 58867 of 1153622 records, EQ = 0.05.

#Exclude problematic records
data_tetrapods_filter_spatialpoints <- data_tetrapods_filter[
  flags$.summary,
  ]

#The flagged records
# dat_fl <- data_teste[!flags$.summary,]

# Remove records that are temporal outliers
flags_temporal <- cf_age(x = data_tetrapods_filter_spatialpoints,
                lon = "decimalLongitude",
                lat = "decimalLatitude",
                taxon = "species", 
                min_age = "year", 
                max_age = "year", 
                value = "flagged")

data_tetrapods_filter_spatial_temp <- data_tetrapods_filter_spatialpoints[
  flags_temporal,
  ]

# Remove records with low coordinate precision
#data_tetrapods_filter_spatial_temp %>% 
#  mutate(Uncertainty = coordinateUncertaintyInMeters / 1000) %>% 
#  ggplot(aes(x = Uncertainty)) + 
#  geom_histogram() +
#  xlab("Coordinate uncertainty in meters") +
#  theme_bw()

# filter coordinates uncertain > 100km
# filter years from 1900
data_tetrapods_filtered <- data_tetrapods_filter_spatial_temp %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters)) %>%
  filter(year >1899)

# Save data
save(data_tetrapods_filtered, # dados limpos
     file = file.path(
       "01_data_cleaned",
       "data_gbif_filtered.RData")
)
end_time <- Sys.time()
print(end_time - start_time) # Time difference of 16.41639 mins
# check data from gbif ----
# check difference per species
length(unique(data_tetrapods_clean$speciesKey)) # 9786 spp
length(unique(data_tetrapods_filtered$speciesKey)) # 7900 spp

# check difference per class
count_species_clean <- data_tetrapods_clean %>%
  distinct(species, .keep_all = TRUE)
table(count_species_clean$class)

count_species_filter <- data_tetrapods_filtered %>%
  distinct(species, .keep_all = TRUE)
table(count_species_filter$class)

# Comparative plots of original vs filtered GBIF datasets
comparative_plot(data_tetrapods_clean,
                 data_tetrapods_filtered,
                 "year",
                 title = "GBIF dataset")

# lat 0
nrow(data_tetrapods_filtered %>%
       filter(decimalLatitude == 0 & decimalLongitude != 0)) 
# long 0
nrow(data_tetrapods_filtered %>%
       filter(decimalLatitude != 0 & decimalLongitude == 0)) 
# 0 long 0 lat
nrow(data_tetrapods_filtered %>%
       filter(decimalLatitude == 0 & decimalLongitude == 0))

nrow(data_tetrapods_filtered %>% 
       filter(is.na(year))) # 0 anos com NA

nrow(data_tetrapods_filtered %>%
       filter(year == 0)) # 0 anos com zero

nrow(data_tetrapods_filtered %>%
       filter(year > 2024)) # 0 registros depois de 2024

nrow(data_tetrapods_filtered %>%
       filter(year < 1899)) # 0 nomes antes de 1900

nrow(data_tetrapods_clean) - nrow(data_tetrapods_filtered) 

data_tetrapods_filtered %>% 
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year))) # without missing data

# BIOTIME FILTER AND CLEAN ----
rm(list=ls()); gc() # clean local enviroment
start_time <- Sys.time()
load(file.path(
  "00_raw_data",
  "biotime_data.RData")
) 

data_biotime_clean <- biotime_data_key_precleaned %>%
  select("speciesKey", "class", "order", "family", "GENUS_SPECIES", 
         "LATITUDE", "LONGITUDE", "YEAR", "STUDY_ID") %>%
  rename(decimalLatitude = LATITUDE, 
         decimalLongitude = LONGITUDE,
         species = GENUS_SPECIES,
         year = YEAR) %>%
  mutate(origin_of_data = "biotime")

names(data_biotime_clean)

data_biotime_clean %>% 
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year))) # without missing data

flags <- clean_coordinates(x = data_biotime_clean, 
                           lon = "decimalLongitude", 
                           lat = "decimalLatitude",
                           species = "species",
                           tests = c("capitals", "centroids", "equal",
                            "gbif","institutions", "outliers",
                            "seas","zeros")) # Flagged 397771 of 1982430 records, EQ = 0.2.

#Exclude problematic records
data_biotime_filter_spatialpoints <- data_biotime_clean[
  flags$.summary,
]

#The flagged records
# dat_fl <- data_teste[!flags$.summary,]
# Remove records that are temporal outliers
flags_temporal <- cf_age(x = data_biotime_filter_spatialpoints,
                         lon = "decimalLongitude",
                         lat = "decimalLatitude",
                         taxon = "species",
                         min_age = "year",
                         max_age = "year",
                         value = "flagged")

data_biotime_filtered <- data_biotime_filter_spatialpoints[
  flags_temporal,
] # TODO arrumar

end_time <- Sys.time()
print(end_time - start_time)

save(data_biotime_filter_spatialpoints,
     file = file.path(
       "01_data_cleaned",
       "data_biotime_filtered.RData")
)
# check data from biotime ----
data_biotime_clean_compare <- data_biotime_clean %>%
  filter(year >1899 & year < 2025)

# quantos registros foram perdidos depois do coordscleaned
nrow(data_biotime_clean_compare) - nrow(data_biotime_filter_spatialpoints)

comparative_plot(data_biotime_clean_compare,
                 data_biotime_filter_spatialpoints,
                 "year", title = "BioTime dataset")

# conferencia
# lat 0
nrow(data_biotime_filter_spatialpoints %>%
       filter(decimalLatitude == 0 & decimalLongitude != 0)) 
# long 0
nrow(data_biotime_filter_spatialpoints %>%
       filter(decimalLatitude != 0 & decimalLongitude == 0)) 
# 0 long 0 lat
nrow(data_biotime_filter_spatialpoints %>%
       filter(decimalLatitude == 0 & decimalLongitude == 0))

nrow(data_biotime_filter_spatialpoints %>% 
       filter(is.na(year))) # 0 anos com NA

nrow(data_biotime_filter_spatialpoints %>%
       filter(year == 0)) # 0 anos com zero

nrow(data_biotime_filter_spatialpoints %>%
       filter(year > 2024)) # 0 registros depois de 2024

nrow(data_biotime_filter_spatialpoints %>%
       filter(year < 1899)) # 0 nomes antes de 1900

# plot time series of biotime
df_filter <- as.data.frame(table(data_biotime_clean$year)) %>%
  rename(Year = Var1) %>%
  mutate(Year = as.numeric(as.character(Year)))

# Criar o gráfico comparativo
ggplot(df_filter, aes(x = Year, y = Freq)) +
  geom_bar(stat = "identity", position = "dodge", width = 2) +
  geom_line(stat = "smooth", linewidth = 1.0, linetype = "dashed") + 
  theme_bw() +
  labs(title = "BioTime", y = "Frequency of occurrences")

# nenhuma coordenada zerada
#nrow(data_biotime_clean %>%
#       filter(decimalLatitude == 0 | decimalLongitude == 0))

# SPLINK FILTER AND CLEAN ----
rm(list=ls()); gc() # clean local enviroment
# start_time <- Sys.time()
start_time <- Sys.time()
load(file.path(
  "00_raw_data",
  "splinks_data.RData")
) 

# names(splink_data_key_precleaned)
data_splink_clean <- splink_data_key_precleaned %>%
  select("speciesKey", "taxonclass", "ordem", "family", "scientificname", 
         "latitude", "longitude", "yearcollected", "coordinateprecision") %>%
  rename(decimalLatitude = latitude, 
         decimalLongitude = longitude,
         class = taxonclass,
         order = ordem,
         species = scientificname,
         year = yearcollected) %>%
  mutate(origin_of_data = "splink") %>%
  mutate(year = na_if(year, "")) %>%
  mutate(coordinateprecision = na_if(coordinateprecision, "")) %>% 
  mutate(year = as.integer(as.character(year))) %>%
  filter(!is.na(year))

#data_tetrapods_filter$countryCode <-  countrycode(
#  data_tetrapods_filter$countryCode,
#  origin =  'iso2c',
#  destination = 'iso3c')

# invalidity coordinates
linhas <- c(19875, 20466, 20880, 21001, 28844, 30254, 33799, 34412, 34414, 34417, 
            34419, 34421, 34423, 34425, 34427, 34428, 34484, 34503, 34505, 34507, 
            34509, 34543, 34544, 50890, 169028, 170394, 170509, 171025, 171026, 
            211771, 211792, 211838, 317962, 344777, 410378, 411561, 411826, 411953, 
            412124, 412125, 413906, 414511, 422472, 424368, 424546, 424563, 426423, 
            426424, 428909, 429430, 429772, 429885, 430079, 445636, 445642, 445728, 
            445729, 445730, 445731, 445732, 445733, 445734, 445735, 445736, 445737, 
            445738, 445739, 445740, 445741, 445742, 445743, 583627)
data_splink_cleaned <-  data_splink_clean[-linhas, ]

flags <- clean_coordinates(x = data_splink_cleaned, 
                           lon = "decimalLongitude", 
                           lat = "decimalLatitude",
                           species = "species",
                           tests = c("capitals", "centroids", "equal",
                            "gbif","institutions", "outliers",
                            "seas","zeros")) # EQ = 0.39.

#Exclude problematic records
data_splink_filter_spatialpoints <- data_splink_cleaned[
  flags$.summary,
]

#The flagged records
# dat_fl <- data_teste[!flags$.summary,]

# Remove records that are temporal outliers
flags_temporal <- cf_age(x = data_splink_filter_spatialpoints,
                         lon = "decimalLongitude",
                         lat = "decimalLatitude",
                         taxon = "species",
                         min_age = "year",
                         max_age = "year",
                         value = "flagged") # Flagged 16743 records

data_splink_filter_spatial_temp <- data_splink_filter_spatialpoints[
  flags_temporal,
]
nrow(data_splink_filter_spatial_temp) # 348817 registros

# Remove records with low coordinate precision
#data_tetrapods_filter_spatial_temp %>% 
#  mutate(Uncertainty = coordinateUncertaintyInMeters / 1000) %>% 
#  ggplot(aes(x = Uncertainty)) + 
#  geom_histogram() +
#  xlab("Coordinate uncertainty in meters") +
#  theme_bw()

# filter coordinates uncertain > 100km
# filter years from 1900
data_splink_filtered <- data_splink_filter_spatial_temp %>%
  filter(coordinateprecision / 1000 <= 100 | is.na(coordinateprecision)) %>%
  filter(year >1899 & year < 2025)
nrow(data_splink_filtered) # 279494

end_time <- Sys.time()
print(end_time - start_time) # 7min57s

save(data_splink_filtered,
  file = file.path(
    "01_data_cleaned",
    "data_splink_filtered.RData")
  )

# check data from splink data ----
# the wrongs years difficult to visualize temporal series tendence 
data_splink_clean_compare <- data_splink_clean %>%
  filter(year >1899 & year < 2025)

comparative_plot(data_splink_clean_compare, data_splink_filtered,
                 "year", title = "splink dataset")

# check difference per species
length(unique(data_splink_filtered$speciesKey)) # 6264 spp
length(unique(data_splink_clean_compare$speciesKey)) # 7252 spp

# Comparative plots of original vs filtered GBIF datasets
comparative_plot(data_splink_clean_compare,
                 data_splink_filtered,
                 "year",
                 title = "splink dataset")
# lat 0
nrow(data_splink_filtered %>%
       filter(decimalLatitude == 0 & decimalLongitude != 0)) 
# long 0
nrow(data_splink_filtered %>%
       filter(decimalLatitude != 0 & decimalLongitude == 0)) 
# 0 long 0 lat
nrow(data_splink_filtered %>%
       filter(decimalLatitude == 0 & decimalLongitude == 0))

nrow(data_splink_filtered %>% 
       filter(is.na(year))) # 0 anos com NA

nrow(data_splink_filtered %>%
       filter(year == 0)) # 0 anos com zero

nrow(data_splink_filtered %>%
       filter(year > 2024)) # 0 registros depois de 2024

nrow(data_splink_filtered %>%
       filter(year < 1899)) # 0 nomes antes de 1900

nrow(data_splink_clean) - nrow(data_splink_filtered) # 318.120

data_splink_filtered %>% 
  summarise(missing_lat = sum(is.na(decimalLatitude)),
            missing_long = sum(is.na(decimalLongitude)), 
            missing_year = sum(is.na(year))) # without missing data

# UNIFYING OCCURENCES DATA ---- 
# TODO aqui vai ser necessario checar quando usarmos dados globais
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

# load data
rm(list=ls()); gc() # clean local enviroment
load(file.path(
  "01_data_cleaned",
  "data_gbif_filtered.RData")
) 
load(file.path(
  "01_data_cleaned",
  "data_biotime_filtered.RData")
)
load(file.path(
  "01_data_cleaned",
  "data_splink_filtered.RData")
)

chaves_perdidas_biotime <- setdiff(
  data_biotime_filter_spatialpoints$speciesKey,
  data_tetrapods_filtered$speciesKey
) 

data_biotime_filter_spatialpoints_sa <- data_biotime_filter_spatialpoints %>% 
  filter(!(speciesKey %in% chaves_perdidas_biotime)) 

chaves_perdidas_splink <- setdiff(
  data_splink_filtered$speciesKey,
  data_tetrapods_filtered$speciesKey
) 

data_splink_filtered_sa <- data_splink_filtered %>% 
  filter(!(speciesKey %in% chaves_perdidas_splink))

glimpse(data_splink_filtered_sa)
glimpse(data_biotime_filter_spatialpoints_sa)

data_occurences_precleaned <- bind_rows(
  data_tetrapods_filtered,
  data_biotime_filter_spatialpoints_sa,
  data_splink_filtered_sa)

save(data_occurences_precleaned, # dados limpos
     file = file.path(
       "01_data_cleaned",
       "data_occurences_cleaned.RData"))

# conference of data ----
table(data_wallacean_knownledge$origin_of_data)

df_count <- data_wallacean_knownledge %>%
  group_by(year, origin_of_data) %>%
  summarise(Occurences = n()) %>%
  ungroup() %>%
  filter(year >1899)

df_total <- df_count %>%
  group_by(year) %>%
  summarise(Total_Occurences = sum(Occurences)) %>%
  filter(year >1899)

ggplot() +
  geom_bar(data = df_count, aes(x = year, y = Occurences, fill = origin_of_data), 
           stat = "identity", position = "stack", width = 0.7, alpha = 0.4) +  # Barras empilhadas
  geom_line(data = df_count, aes(x = year, y = Occurences, color = origin_of_data, group = origin_of_data), 
            size = 1.0, linetype = "dashed", stat = "smooth") +  # Linhas para cada origem de dados
  #geom_smooth(data = df_total, aes(x = year, y = Total_Occurences), 
  #            method = "loess", se = FALSE, color = "black", size = 1.0, linetype = "dashed") +
  theme_bw() +
  labs(title = "Série Temporal de ocorrências por origem dos dados", 
       x = "Ano", 
       y = "Número de Ocorrências",
       fill = "Origem dos Dados",
       color = "Origem dos Dados") +  # Sincroniza legendas
  scale_fill_brewer(palette = "Set1") +  # Paleta para barras
  scale_color_brewer(palette = "Set1")

nrow(data_wallacean_knownledge) # 2226686 total occurences 
table(data_wallacean_knownledge$class) # occurences per classes 

# REMOVE DUPLICATES AND FILTER POINTS BY RANGE POLYGONS ----
rm(list=ls()); gc() # clean local enviroment
start_time <- Sys.time()
# shapefile
load(file = file.path(
  "00_raw_data",
  "shapefiles_data.RData"))
# occurences data
load(file = file.path(
  "01_data_cleaned",
  "data_occurences_cleaned.RData"))

# plot points without polygons filter
#ggplot() +
#  coord_fixed() +
#  borders("world", colour = "gray50", fill = "gray50") +
#  geom_point(data = data_occurences_precleaned,
#             aes(x = decimalLongitude, y = decimalLatitude),
#             colour = "darkred",
#             size = 0.5) +
#  theme_bw()

# Remove duplicates
duplicated_flags <- cc_dupl(
  data_occurences_precleaned,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  species = "speciesKey",
  value = "flagged"
)

data_occurences_precleaned_duplicate <- data_occurences_precleaned[
  duplicated_flags == TRUE, ]

# deixar apenas as especies que tem poligonos
data_wallacean_knownledge_sa <- data_occurences_precleaned_duplicate %>%
  filter(speciesKey %in% tetrapod_shapefile$speciesKey)
# deixar apenas os poligonos que tem dados de distribuicao
data_tetrapods_sa <- tetrapod_shapefile %>% 
  filter(speciesKey %in% data_occurences_precleaned_duplicate$speciesKey)
# tem que ter o mesmo numero de especies (speciesKey)
length(unique(data_wallacean_knownledge_sa$speciesKey)) # 6548 spp
length(unique(data_tetrapods_sa$speciesKey)) # 6548 spp

# run cc_iucn()
# terra::vect(data_tetrapods_sa) se funcionar vai dar certo
# Tentar converter todas as geometrias para MULTIPOLYGON
data_tetrapods_sa_transf <- st_cast(data_tetrapods_sa,
 "MULTIPOLYGON")
# Verificar se a conversão funcionou
geometry_types <- st_geometry_type(data_tetrapods_sa_transf)
table(geometry_types) 

# check if list_occurences_clean is empty before run
unique_species_keys <- unique(data_wallacean_knownledge_sa$speciesKey)
list_occurences_clean <- data.frame()
species_keys_with_errors <- c()

for (species_key in unique_species_keys) {
  # Usar tryCatch para capturar erros
  tryCatch({
    # Selecionar o polígono correspondente à espécie
    especie_poligono <- data_tetrapods_sa_transf %>%
      filter(speciesKey == species_key)
    
    subset_data <- data_wallacean_knownledge_sa %>%
      filter(speciesKey == species_key)
    
    # Aplicar a função cc_iucn
    range_flags <- cc_iucn(
      x = subset_data,
      range = especie_poligono,
      species = "speciesKey",
      lon = "decimalLongitude",
      lat = "decimalLatitude",
      value = "flagged"
    )
    
    # Filtrar os dados usando os valores retornados por range_flags
    subset_data_filtered <- subset_data[range_flags == TRUE, ]
    
    # Adicionar os dados filtrados à lista final
    list_occurences_clean <- rbind(list_occurences_clean, subset_data_filtered)
    
  }, error = function(e) {
    # Se ocorrer um erro, salvar o species_key no vetor
    species_keys_with_errors <<- c(species_keys_with_errors, species_key)
    # Opcional: Printar o erro ou a species_key para depuração
    message(paste("Erro no speciesKey:", species_key, "-", e$message))
  })
}
end_time <- Sys.time()
print(end_time - start_time) # Time difference of 1:30+-

# Verificar speciesKey que tiveram erros
print("SpeciesKey com erro:")
print(species_keys_with_errors)

#data_tetrapods_sa_transf <- data_tetrapods_sa %>%
#    filter(speciesKey %in% species_keys_with_errors)

nrow(list_occurences_clean) # 267.722
nrow(data_wallacean_knownledge_sa) # 109931 excluidos pelo poligono
length(unique(list_occurences_clean$speciesKey)) # 5726 spp
length(unique(data_wallacean_knownledge_sa$speciesKey)) # 6548 spp
table(list_occurences_clean$class)
# Compare
plot_sujo <- ggplot() +
  coord_fixed() +
  borders("world", colour = "gray50", fill = "gray50") +
  geom_point(data = data_wallacean_knownledge_sa,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred",
             size = 0.5) +
  theme_bw()

plot_clean <- ggplot() +
  coord_fixed() +
  borders("world", colour = "gray50", fill = "gray50") +
  geom_point(data = list_occurences_clean,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred",
             size = 0.5) +
  theme_bw()

# Save 
save(list_occurences_clean,
     species_keys_with_errors,
     file = here(
       "01_data_cleaned",
       "list_occurences_clean.RData")
)

# NESTED DATAFRAME ----
rm(list=ls()); gc() # clean local enviroment
load(file = here(
       "01_data_cleaned",
       "list_occurences_clean.RData")
)
load(file.path(
  "00_raw_data",
  "trait_data.RData")
) # match taxonomic and trait information with speciesKey

tetrapods_key_species <- trait_data %>%
  select("speciesKey","scientificName","Class", "Order","Family", 
  "YearOfDescription","BodyLength_mm","ImputedLength","BodyMass_g",
  "ImputedMass","Diu","Noc","Nocturnality", "Fos", "Ter", "Aqu", "Arb", "Aer",
  "ImputedHabitat","MajorHabitatSum","ImputedMajorHabitat","RangeSize",
   "HumanDensity", "AssessedStatus") %>%
  distinct(speciesKey, .keep_all = TRUE) # 117 speciesKey duplicated
anyDuplicated(tetrapods_key_species$speciesKey) # ok

data_tetrapods_nested <- list_occurences_clean %>%
  group_by(speciesKey) %>%  # Agrupa pelo speciesKey
  nest() %>%
  rename(event_table = data)
anyDuplicated(data_tetrapods_nested$speciesKey) # ok

data_wallacean_nested <- left_join(
  data_tetrapods_nested,
  tetrapods_key_species,
  by = "speciesKey"
) %>% 
  mutate(count_events = NA) %>%
  relocate(event_table,count_events, .after = Family) %>%
  arrange(Class, Order, Family)

head(data_wallacean_nested)

# check data
table(data_wallacean_nested$Class)
# Amphibia     Aves   Mammalia Reptilia
#     789      2765      893     1279
table(data_wallacean_nested$Order)

anyDuplicated(data_wallacean_nested$scientificName) # nenhum nome repetido
anyDuplicated(data_wallacean_nested$speciesKey)

for (i in 1:nrow(data_wallacean_nested)) {
  data_wallacean_nested[i,"count_events"] <- nrow(
    data_wallacean_nested$event_table[[i]])
  
  data_wallacean_nested$event_table[[i]] <- 
    data_wallacean_nested$event_table[[i]][
      order(data_wallacean_nested$event_table[[i]]$year), 
    ]
}

save(data_wallacean_nested, 
     file = here(
       "01_data_cleaned",
       "data_wallacean_nested.RData")
)

# Explorar os dados
ggplot(data_wallacean_nested, aes(x = log(RangeSize), y = log(count_events))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Adiciona uma linha de tendência
  labs(x = "Comprimento do Corpo (mm)", y = "Número de Registros (count_events)", title = "Relação entre Comprimento do Corpo e Registros") +
  theme_minimal()

ggplot(data_wallacean_nested, aes(x = AssessedStatus, y = log(count_events))) +
  geom_boxplot(fill = "lightblue", colour = "darkblue") +
  theme_minimal() +
  labs(title = "Comparação de quantidade de registros por categoria da IUCN",
       x = "Categoria da IUCN",
       y = "Quantidade de registros (event_counts)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
